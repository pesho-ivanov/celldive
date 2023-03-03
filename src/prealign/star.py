#!/usr/bin/python3
import argparse
import configparser 
import subprocess
import argparse
import argcomplete

import tools
from tools import VerifyFile
from tools import create_dir
from tools import match_dir

# TODO: http://stackoverflow.com/questions/3609852/which-is-the-best-way-to-allow-configuration-options-be-overridden-at-the-comman
def get_args():
    parser = argparse.ArgumentParser()
    argcomplete.autocomplete(parser)

    subparsers = parser.add_subparsers(dest='mode')
    subparsers.required = True

    # FOR ALL
    parser.add_argument("-c", "--config",
                        help="Config file for bioinf libs and data",
                        action=VerifyFile(), default=tools.DEFAULT_CONF)
    parser.add_argument('-t', '--threads',
                        help="Number of threads", type=int, default=1)
    parser.add_argument('-o', '--output_dir',
                        help="Output dir", action=VerifyFile(), type=create_dir)
  # parser.add_argument('fn',action=VerifyExt({'csv','txt'}))
#    parser.add_argument('-s', '--sample_percent',
#            help="Sample this much from every fastq file.", type=float)

    # INDEX BUILDING
    parser_build = subparsers.add_parser(
        'build', help='Build an index based on genome files and annotation.')
    parser_build.add_argument(
        '--input_genome', #type=lambda fn: match_fasta(fn, gz=False),
        help="Input genome/transcriptome fasta file.", required=True)
    parser_build.add_argument(
        '--input_gtf', #type=lambda fn: tools.match_file_type(fn, '.gtf'),
        help="Annotations file for the genome/transcriptome.")
#    parser_build.add_argument('-o', '--output_dir',
#            type=create_dir,
#            help="Directory to be created for the STAR index.")
    parser_build.add_argument(
        '--overhang',
        help="Overhang parameter to be passed to STAR.", type=int, default=100)

    # SINGLE RUN
    parser_single = subparsers.add_parser(
        'singlerun', help='Run over a single input.')
    parser_single.add_argument(
        '-1', '--fastq1', help='One fastq(.gz) file',
        action=VerifyFile({'fq', 'fastq'}), required=True)
    parser_single.add_argument(
        '-2', '--fastq2',
        help='The pairing fastq(.gz) file',
        action=VerifyFile({'fq', 'fastq'}), required=False)

    # BATCH RUN
    parser_batch = subparsers.add_parser(
        'batchrun', help='Run over a collection of inputs.')
    parser_batch.add_argument(
        'inputdir', help='A directory with fastq(.gz) files',
        action=VerifyFile(), type=match_dir)

    args = parser.parse_args()
    return args

def build_index(input_genome, input_gtf, output_dir, overhang, threads,
        STAR_BIN):
    cmd = [STAR_BIN,
        '--runMode',            'genomeGenerate',
        '--genomeFastaFiles',   input_genome,
        '--genomeDir',          output_dir,
        '--runThreadN',         str(threads),
    ]
    if input_gtf:
        cmd.extend(['--sjdbGTFfile', input_gtf])
        cmd.extend(['--sjdbOverhang', str(overhang)])
    if is_gz(input_genome):
        cmd.extend(['--readFilesCommand', 'zcat'])
    tools.run_cmd(cmd)

def run_star(fastqs, output_dir, threads,
        STAR_BIN, STAR_INDICES_DIR, GTF_FILE):
    assert len(fastqs) <= 2
    cmd = [STAR_BIN,
        '--genomeDir',          STAR_INDICES_DIR,
        '--sjdbGTFfile',        GTF_FILE,
        '--readFilesIn',        ' '.join(fastqs),
        '--outFileNamePrefix',  output_dir,
        '--runThreadN',         str(threads),
    ]
    if are_gz(fastqs):
        cmd.extend(['--readFilesCommand', 'zcat'])
    tools.call_cmd(cmd)

def batchrun_star():
    return

if __name__ == "__main__":
    args = get_args()
    LIBS, DATA = tools.get_config(args.config)

    if args.mode == 'build':
        build_index(args.input_genome, args.input_gtf, args.output_dir,
                    args.overhang, args.threads,
                    LIBS['star_bin'])
    elif args.mode == 'singlerun':
        fastqs = [args.fastq1]
        if args.fastq2:
            fastqs.append(args.fastq2)
        run_star(fastqs, args.output_dir, args.threads,
                 LIBS['star_bin'], DATA['star_indices_dir'], DATA['gtf_file'])
    elif args.mode == 'batchrun':
        batchrun_star()
