#!/bin/bash
# Runs bowtie2 on all the reads

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
CONFIG_FILE=/home/pesho/libs/tracer-0.4.0/tracer-human.conf
CORES=25

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir $OUTPUT_DIR

echo "### run-human.sh ###"
echo "Processing all .fastq(.gz) files in ${INPUT_DIR}"
echo "Output in ${OUTPUT_DIR}"

COUNTER=0
# Run tracer assemble for every cell
for longfn in $INPUT_DIR/*.fastq*
do
   fn=${longfn##*/}
   base=${fn%.fastq.gz}
   if [ ${base:(-2)} == "_1" ]; then
     echo "`date`"
     echo $base
     ~/libs/kallisto-0.43.0/kallisto quant -i /bigcode/pesho/resources/transcriptomes/Homo_sapiens.GRCh38.rel79.cdna.all.idx -o $OUTPUT_DIR/$base $longfn ${longfn/_1.fastq.gz/_2.fastq.gz} &
     let COUNTER=COUNTER+1
     if (( $COUNTER % $CORES == 0 )); then wait; fi
   else 
     echo 'WARNING: ' $base ' not processed!'
   fi
done

# waiting for all parallel cells
wait

