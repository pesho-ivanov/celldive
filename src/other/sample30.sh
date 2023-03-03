#!/bin/bash
# Samples 30 random cells

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
SAMPLE_SIZE=30

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
ls $INPUT_DIR > sample_from.txt
shuf -n $SAMPLE_SIZE sample_from.txt > sample.txt

for f in `cat sample.txt`; do
  #if [ $f ne 'filtered_TCR_summary' ]; do
  echo $f
  cp -r $INPUT_DIR/$f .
  #fi
done
