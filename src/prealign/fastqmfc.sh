#!/bin/bash

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
CORES=25

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir -p $OUTPUT_DIR

echo "### Trimming ###"
echo "Processing all .fastq(.gz) files in ${INPUT_DIR}"
echo "Output in ${OUTPUT_DIR}"

COUNTER=0
# Run tracer assemble for every cell
for longfn in $INPUT_DIR/*.fastq*
do
    fn=${longfn##*/}
    firstbase=${fn%.fastq.gz}
    if [ ${firstbase:(-1)} == "1" ];
    then
        base=${firstbase::(-1)}
        second=${longfn/1.fastq.gz/2.fastq.gz}
        echo $longfn $second

        out1=$OUTPUT_DIR/${base}_1.fastq.gz
        out2=$OUTPUT_DIR/${base}_2.fastq.gz

        # Usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...]
        fastq-mcf \
            ~/libs/trimmomatic-0.36/adapters/NexteraPE-PE.fa \
            $longfn $second -o $out1 -o $out2   &   # DONT WAIT

        let COUNTER=COUNTER+1
        if (( $COUNTER % $CORES == 0 )); then wait; fi
   fi
done
wait

fastqc -t $CORES $OUTPUT_DIR/*
multiqc $OUTPUT_DIR/ -o $OUTPUT_DIR -p
