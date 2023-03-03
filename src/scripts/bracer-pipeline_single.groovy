// Exec: bpipe run -d out/ pipeline.groovy in/*.fastq.gz
// Sample pipeline: http://bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html
//     TopHat (on GRCh38 and RefSeq human gene annotation) for alignment
//   + featureCounts, cufflinks for counts and FPKM
//   + RSEM (to transcriptome) for TPM
//   + FastQC and Picard for QC

// REF="~/hg19/gatk.ucsc.hg19.fasta"

// NOTE: for single ended data
// is WORK IN PROGRESS

// current settings for pisch
// note: activate environment for bracer with
//       conda activate bracer
KALLISTO="/home/pesho/libs/kallisto-0.43.1/kallisto"
//index=/data/bio/depend/transcriptomes/human/Homo_sapiens.GRCh38.rel79.cdna.all.idx
//index=/data/bio/depend/transcriptomes/qc/ArrayControl.idx
INDEX="/data/bio/depend/transcriptomes/human/Homo_sapiens+ArrayControl.GRCh38.rel79.cdna.all.idx"
MERGE_TSVS_SCRIPT="/home/pisch/CellDive/src/scripts/merge-tsv-to-table.sh"
BRACER_BIN="/home/pisch/anaconda3/envs/bracer/bin/bracer"
CONFIG_FILE="/home/pisch/CellDive/src/scripts/bracer-human.conf"

about title: "Bracer"

bracer_run = {
	output.dir = 'bracer'

        //produce("quantify_kallisto.log", "abundance.tsv") {
          from("fastq*") {
            exec """
            $BRACER_BIN assemble \
                -c $CONFIG_FILE -s Hsap -p 1 --single_end \
                $branch $output.dir $input1
            """
          }
        //}
}

bracer_summarize = {
	output.dir = 'bracer'
        exec "$BRACER_BIN summarise -c $CONFIG_FILE $output.dir"
}

run {
	"%_*.fastq*" * [ bracer_run ] \
        + bracer_summarize 
}
