#!/bin/bash
# Run Kallisto

#echo ${##[@]}
if [ $# -lt 3 ]; then
  echo 'Please give two arguments:'
  echo '  -- output counts spreadsheet file'
  echo '  -- output counts spreadsheet file'
  echo '  -- input dir from Kallisto'
  exit 1
fi

OUTPUT_FILE_COUNTS=$1;
OUTPUT_FILE_TPM=$2;  
#shift 2
INPUTS=("$@")
INPUTS=("${INPUTS[@]:2}")
 
echo 'inputs: ' $OUTPUT_FILE_COUNTS $OUTPUT_FILE_TPM
if [ ${#INPUTS[@]} -lt 1 ]; then
  echo 'No input files'
  return 1
fi

#if [ -z $OUTPUT_FILE_TMP ]; then
#  OUTPUT_FILE_TPM=$INPUT_DIR/expressions_tpm.tsv
#fi

#if [ -z $OUTPUT_FILE_COUNTS ]; then
#  OUTPUT_FILE_COUNTS=$INPUT_DIR/expressions_counts.tsv
#fi

#echo ${#INPUTS[@]}
#echo 'all input files:  ' ${INPUTS}
#echo 'first input file: ' ${INPUTS[0]}

cat ${INPUTS[0]} | awk '{print $1}' > $OUTPUT_FILE_TPM
cat ${INPUTS[0]} | awk '{print $1}' > $OUTPUT_FILE_COUNTS

for file in ${INPUTS[@]}
do
  echo 'merging to the table: ' $file
  CELL=`echo $file | sed 's/.*kallisto\/\(.*\)\/.*/\1/'`

  echo 'cell: ' $CELL

#  cat $OUTPUT_FILE_COUNTS | tail
  paste $OUTPUT_FILE_COUNTS <(cat $file | awk '{print $4}' | sed -e "1s/.*/$CELL/") > tmp
  mv tmp $OUTPUT_FILE_COUNTS
  paste $OUTPUT_FILE_TPM <(cat $file | awk '{print $5}' | sed -e "1s/.*/$CELL/") > tmp
  mv tmp $OUTPUT_FILE_TPM
done

echo 'output files: ' $OUTPUT_FILE_COUNTS $OUTPUT_FILE_TPM
