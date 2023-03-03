// Exec: bpipe run -d out/ pipeline.groovy in/*.fastq.gz
// Sample pipeline: http://bioconductor.org/packages/release/data/experiment/vignettes/scRNAseq/inst/doc/scRNAseq.html
//     TopHat (on GRCh38 and RefSeq human gene annotation) for alignment
//   + featureCounts, cufflinks for counts and FPKM
//   + RSEM (to transcriptome) for TPM
//   + FastQC and Picard for QC

// dependency check: https://projects.nbis.se/projects/bpipe/wiki/Bpipe_developer_notes

// REF="~/hg19/gatk.ucsc.hg19.fasta"


//--- settings for Pesho
//KALLISTO="/home/pesho/libs/kallisto-0.43.1/kallisto"
//index=/data/bio/depend/transcriptomes/human/Homo_sapiens.GRCh38.rel79.cdna.all.idx
//index=/data/bio/depend/transcriptomes/qc/ArrayControl.idx
INDEX="/data/bio/depend/transcriptomes/human/Homo_sapiens+ArrayControl.GRCh38.rel79.cdna.all.idx"
MERGE_TSVS_SCRIPT="/home/pesho/repos/bioviz/src/scripts/merge-tsv-to-table.sh"

MULTIQC="multiqc"
KALLISTO_OUT_DIR = "kallisto"
KALLISTO="kallisto"


//--- settings for pisch
// note: needs to run in the conda environment 'multiqc'
// note2: currently, bpipe needs to be launched from the output target folder
//KALLISTO="/home/pisch/.local/bin/kallisto"
// Hs
//INDEX="/data/bio/depend/transcriptomes/human/Homo_sapiens+ArrayControl.GRCh38.rel79.cdna.all.idx"
// Mm
//INDEX="/data/bio/depend/transcriptomes/other/Mus_musculus.GRCm38.rel79.cdna.all.idx"
//MERGE_TSVS_SCRIPT="/home/pisch/CellDive/src/scripts/merge-tsv-to-table.sh"



about title: "Aligning, quantification"


quantify_kallisto = {
	output.dir = KALLISTO_OUT_DIR + "/" + branch

        produce("quantify_kallisto.log", "abundance.tsv") {
        from("*", "*") {
                exec "$KALLISTO quant -i $INDEX -o ${output.dir} -t $threads --bias $input1 $input2 2> $output"
          }
        }
}

merge_kallistos_into_table = {
        output.dir = KALLISTO_OUT_DIR

        produce("kallisto_counts.tsv", "kallisto_tpm.tsv") {
                exec "$MERGE_TSVS_SCRIPT $output1 $output2 $inputs"
        }
        //exec "$MERGE_TSVS_SCRIPT kallisto_counts.tsv kallisto_tmp.tsv $inputs"
}

multiqc = {
        output.dir = "multiqc"
        //multiqc $OUTPUT_DIR/ -o $OUTPUT_DIR -p

        exec "$MULTIQC -p $KALLISTO_OUT_DIR -o $output.dir"
}

run {
	"%_%_%_%_*" * [ quantify_kallisto ] \
        + "*abundance.tsv" * [ merge_kallistos_into_table ] \
        + multiqc
}
