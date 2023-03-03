#!/bin/bash
# Runs kallisto on all the reads
#
# run in the multiqc python environment

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- dir with fastq.gz files to look into (Texts)'
  echo '  -- output dir for the .sam files'
  exit 1
fi

kallisto=/home/pisch/.local/bin/kallisto
index=/data/bio/depend/transcriptomes/other/Mus_musculus.GRCm38.rel79.cdna.all.idx
merge_expr=/home/pisch/CellDive/src/scripts/merge-tsvs.sh

FILES=$1
OUTPUT_DIR=$2
CORES=25

if [ ! -d $OUTPUT_DIR ]; then
    mkdir -p $OUTPUT_DIR
fi

if [ ! -d $FILES ]; then
    echo 'No such dir' $FILES
    exit 1
fi

COUNTER=0
for first in $FILES/*.fastq.gz; do
   fn=${first##*/}
   base=${fn%.fastq.gz}
   echo "Processing $first $fn $base"

   if [ ${base:(-1)} == "1" ]; then
     bb=${base::(-1)}
     second=${first/1.fastq.gz/2.fastq.gz}

     [ -f ${first} ] || echo "File $first not found"
     [ -f ${second} ] || echo "File $second not found"

     #$kallisto quant --threads=30 -i $index -o $OUTPUT_DIR/$bb --bias $first $second
     $kallisto quant -i $index -o $OUTPUT_DIR/$bb --bias $first $second 2> $OUTPUT_DIR/kallisto.log &
     let COUNTER=COUNTER+1
     if (( $COUNTER % $CORES == 0 )); then wait; fi
   else
     if [ ! ${base:(-1)} == "2" ]; then
       echo "Strange filename: " $first
     fi
   fi
done
wait

$merge_expr $OUTPUT_DIR $OUTPUT_DIR/expressions_counts.tsv $OUTPUT_DIR/expressions_tpm.tsv
multiqc $OUTPUT_DIR/ -o $OUTPUT_DIR -p
