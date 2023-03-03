#!/bin/bash
# Run Kallisto

if [ $# -ne 3 ]; then
  echo 'Please give three arguments:'
  echo '  -- input dir from Kallisto'
  echo '  -- output counts spreadsheet file'
  echo '  -- output TPM spreadsheet file'
  exit 1
fi
INPUT_DIR=$1
OUTPUT_FILE_COUNTS=$2  #$INPUT_DIR/expressions_counts.tsv
OUTPUT_FILE_TPM=$3  #$INPUT_DIR/expressions_tpm.tsv

#if [ -z $OUTPUT_FILE_TMP ]; then
#  OUTPUT_FILE_TPM=$INPUT_DIR/expressions_tpm.tsv
#fi

#if [ -z $OUTPUT_FILE_COUNTS ]; then
#  OUTPUT_FILE_COUNTS=$INPUT_DIR/expressions_counts.tsv
#fi

FIRST_DIR=`ls -d $INPUT_DIR/* | sort -n | head -1`
cat $FIRST_DIR/abundance.tsv | awk '{print $1}' > $OUTPUT_FILE_TPM
cat $FIRST_DIR/abundance.tsv | awk '{print $1}' > $OUTPUT_FILE_COUNTS

for file in $INPUT_DIR/*/*.tsv
do
  echo 'merging to the table: ' $file
  #batch3: CELL=`echo $file | sed 's/.*-C1_\([0-9]*\).*/\1/'`
  #batch1:
  #CELL=`echo $file | sed 's/.*_\([0-9]*\)_.*/\1/'`
  CELL=`echo $file | sed 's/[^0-9]*\([0-9]*\)_.*/\1/'`
  #batch2: CELL=`echo $file | sed 's/.*_\([0-9C]*\)\/abundance.tsv/\1/'`
  #CELL=`echo $file | sed 's/\.\/\([0-9C]*\)\/abundance.tsv/\1/'`

  paste $OUTPUT_FILE_COUNTS <(cat $file | awk '{print $4}' | sed -e "1s/.*/$CELL/") > tmp
  mv tmp $OUTPUT_FILE_COUNTS
  paste $OUTPUT_FILE_TPM <(cat $file | awk '{print $5}' | sed -e "1s/.*/$CELL/") > tmp
  mv tmp $OUTPUT_FILE_TPM
done

echo 'output files: ' $OUTPUT_FILE_COUNTS $OUTPUT_FILE_TPM
