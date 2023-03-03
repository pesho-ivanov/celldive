// Usage:
//   1) go to the newly created output directory out_dir
//   2) execute
//         bpipe run -n 28 <PATH_TO_SCRIPT>/trim-pipeline.groovy <PATH_TO_SAMPLE>/*
//   3) assume that the sample files are <UNIQUE_NAME>_R(1|2).fastq*

FASTQC="fastqc"
MULTIQC="multiqc"
//ADAPTERS_FILE="~/libs/trimmomatic-0.36/adapters/NexteraPE-PE.fa"
TRIMMED_DIR="trimmed"
POLYT="TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
TRIM_GALORE="trim_galore"
FASTQC_RAW_DIR="fastqc-raw"
FASTQC_TRIMMED_DIR="fastqc-trimmed"
MULTIQC_DIR="multiqc"

about title: "Trimming"

trim_galore_pair_end = {
        output.dir = TRIMMED_DIR

	from("fastq*", "fastq*") {
		exec "$TRIM_GALORE --nextera --paired $input1 $input2 -o ${output.dir}"
//		exec "$TRIM_GALORE -a "$POLYT" --paired $input1 $input2"
	}
	forward input1, input2
}

fastqc_over_pair = {
        output.dir = FASTQC_RAW_DIR
	exec "$FASTQC -t $threads $input1 $input2 -o ${output.dir}"
	forward input1, input2
}

fastqc_over_trimmed = {
        output.dir = FASTQC_TRIMMED_DIR
	exec "$FASTQC -t $threads $TRIMMED_DIR/*.gz -o ${output.dir}"
	forward output.dir
}

multiqc = {
        output.dir = MULTIQC_DIR
        exec "$MULTIQC $input --title trimmed --comment Trimmed -o $output.dir"
}

run {
	"%_R*.fastq*" * [ trim_galore_pair_end + fastqc_over_pair ] + fastqc_over_trimmed + multiqc
}
