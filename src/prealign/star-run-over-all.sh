#!/bin/bash

if [ $# -ne 1 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- genome dir'
  echo '  -- output dir'
  exit 1
fi

GENOME_DIR=$1

for longfn in $INPUT_DIR/*.fastq*
do
    fn=${longfn##*/}
    firstbase=${fn%.fastq.gz}
    if [ ${firstbase:(-2)} == "_1" ];
    then
        base=${firstbase::(-2)}
        second=${longfn/_1.fastq.gz/_2.fastq.gz}

STAR --runThreadN 26 --runMode ? \
  --genomeDir $GENOME_DIR \
  --readFilesIn  \
  --genomeFastaFiles $GENOME_DIR/sequence/GRCh38_r77.all.fa \
  --sjdbGTFfile $GENOME_DIR/annotation/Homo_sapiens.GRCh38.77.gtf \
  --sjdbOverhang 100
