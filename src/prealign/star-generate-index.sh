#!/bin/bash

if [ $# -ne 1 ]; then
  echo 'Please give two arguments:'
  echo '  -- genome dir'
  exit 1
fi

GENOME_DIR=$1
mkdir -p $GENOME_DIR/star_indices_overhang100
STAR --runThreadN 26 --runMode genomeGenerate \
  --genomeDir $GENOME_DIR/star_indices_overhang100/ \
  --genomeFastaFiles $GENOME_DIR/sequence/GRCh38_r77.all.fa \
  --sjdbGTFfile $GENOME_DIR/annotation/Homo_sapiens.GRCh38.77.gtf \
  --sjdbOverhang 100
