#!/bin/bash
# Runs bowtie2 on all the reads

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- dir with fastq.gz files to look into (Texts)'
#  echo '  -- path to bowtie index (to 6 files with diff extensions)'
  echo '  -- output dir for the .sam files'
  exit 1
fi

FILES=$1
PATTERNS=/data/bio/depend/genomes/hg19
OUTPUT_DIR=$2

for longfn in $FILES/*.fastq.gz
do
   fn=${longfn##*/}
   base=${fn%.fastq.gz}
   echo "Processing $base"

   /home/pesho/libs/bowtie-2.2.9/bowtie2 \
    --no-unal -p 20 -k 1 --np 0 -t \
    --rdg 1,1 --rfg 1,1 -x $PATTERNS \
    -U ${longfn} \
    --al $OUTPUT_DIR/${base}.fasta \
    -S $OUTPUT_DIR/${base}.sam 2>&1 | tee $OUTPUT_DIR/$base.log
done
