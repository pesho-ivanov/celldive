// Exec: bpipe run -n <n_cores> --dir out/ --param sample_id=<sample_id> pipeline.groovy in/*.fastq.gz
//
// Sample pipeline: http://bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html
//     TopHat (on GRCh38 and RefSeq human gene annotation) for alignment
//   + featureCounts, cufflinks for counts and FPKM
//   + RSEM (to transcriptome) for TPM
//   + FastQC and Picard for QC

// bpipe run -n 5 -p reference=./reference.fa -p RMDUP="yes" -p REALIGN="yes" -p BAQ="yes" -p SPLIT="%_*.fastq" -p ALNMETHOD="mem" -p SUFFIX=".fastq" ./bwaWrapper.groovy ./*_1.fastq ./*_2.fastq
// run { check_config + bwa_index + "${SPLIT}" * [bwa] + "%.sam" * [sam_view + sam_flagstat, sam_view_unmap, sam_view_notpp, pp] + last

// REF="~/hg19/gatk.ucsc.hg19.fasta"

//KALLISTO="/home/pesho/libs/kallisto-0.43.1/kallisto"
//index=/data/bio/depend/transcriptomes/human/Homo_sapiens.GRCh38.rel79.cdna.all.idx
//index=/data/bio/depend/transcriptomes/qc/ArrayControl.idx
//INDEX="/data/bio/depend/transcriptomes/human/Homo_sapiens+ArrayControl.GRCh38.rel79.cdna.all.idx"
//MERGE_TSVS_SCRIPT="/home/pesho/repos/bioviz/src/scripts/merge-tsv-to-table.sh"

// for pesho
//TRACER_BIN="/home/pesho/libs/tracer-0.5.1/tracer"
//TRACER2JSON="/home/pesho/repos/bioviz/src/data2json/tracer2json.py"
//CONFIG_FILE="/home/pesho/repos/bioviz/src/scripts/tracer-human.conf"

// for pisch
// note: needs to run in activated tracer environment
TRACER_BIN="/home/pisch/anaconda3/envs/tracer/bin/tracer"
TRACER2JSON="/home/pisch/CellDive/src/data2json/tracer2json.py"
CONFIG_FILE="/home/pisch/CellDive/src/scripts/tracer-human_pisch.conf"

about title: "Tracer"
//SUBDIR = "SUBDIR"

tracer_run = {
	//output.dir = "${OUTDIR}/tracer"
        output.dir = "tracer"

        //produce("filtered_TCR_seqs/filtered_TCRs.txt", "abundance.tsv") {
          from("fastq*", "fastq*") {
            exec """
            $TRACER_BIN assemble \
                -c $CONFIG_FILE -s Hsap -p 1 \
                $input1 $input2 $branch $output.dir
            """
          }
        //}
}

tracer_summarize = {
	//output.dir = "${OUTDIR}/tracer"
        output.dir = "tracer"
	exec "$TRACER_BIN summarise --no_network -c $CONFIG_FILE $output.dir"
}

tracer_to_json = {
	output.dir = "tracer"
	exec "python $TRACER2JSON $output.dir $sample_id"
}

run {
	// typical pattern
	//"plate*_%_R*.fastq*" * [ tracer_run ] \
        //+ tracer_summarize \
	//+ tracer_to_json
	
	// special pattern for some batches (first batches)
	//"%_*.fastq*" * [ tracer_run ] \
	//+ tracer_summarize \
	//+ tracer_to_json 

	// new batches
	"*.A-Plate_*_%_R*.fastq*" * [ tracer_run ] \
	+ tracer_summarize \
	+ tracer_to_json
}
