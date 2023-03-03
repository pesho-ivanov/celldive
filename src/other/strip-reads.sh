#!/bin/bash
# Runs bowtie2 on all the reads

if [ $# -ne 3 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  echo '  -- new reads length (leaving the prefix of every read)'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
NEW_LEN=$3

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir -p $OUTPUT_DIR

echo "### `basename "$0"`  ###"
echo "Stripping all reads in .fastq(.gz) files in ${INPUT_DIR} to the first ${NEW_LEN}bp"
echo "Output written in ${OUTPUT_DIR}"

# Run tracer assemble for every cell
for longfn in $INPUT_DIR/*.fastq
do
   fn=${longfn##*/}
   base=${fn%.fastq}
   echo "${fn}..."
   cat $longfn | cut -c -$NEW_LEN > $OUTPUT_DIR/$fn
done

for longfn in $INPUT_DIR/*.fastq.gz
do
   fn=${longfn##*/}
   base=${fn%.fastq.gz}
   echo "${fn}..."
   zcat $longfn | cut -c -$NEW_LEN > $OUTPUT_DIR/${fn%.gz}
done
