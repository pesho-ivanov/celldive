#!/bin/bash
# Merges the .fasta.gz files for each cell

if [ $# -ne 2 ]; then
  echo 'Please give two arguments: input and output directories'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
EXT='fastq'
LOG_FILE=$OUTPUT_DIR/merge-fsq.log

if [ ! -r $INPUT_DIR ]; then
  echo "Failure: The input dir $INPUT_DIR does not exist or is not readable"
  exit 1
fi

if [ -d $OUTPUT_DIR ]; then
  if [ ! -w $OUTPUT_DIR ]; then
    echo "Failure: The output dir $OUTPUT_DIR is not writable"
    exit 1
  fi
  if [ "$(ls -A $OUTPUT_DIR)" ]; then
    echo "Failure: The output dir $OUTPUT_DIR is not empty"
    exit 1
  fi
else
  mkdir -p $OUTPUT_DIR
fi

echo "Merging multiplexed files from $INPUT_DIR into $OUTPUT_DIR..." | tee $LOG_FILE

for longfn in $INPUT_DIR/*.$EXT*
do
   fn=${longfn##*/}
   base=${fn%.$EXT*}

   if [[ $base =~ ^([^_]+_[^_]+_[^_]+) ]]; 
   then 
     merged_file=${BASH_REMATCH[1]}.$EXT
   else 
     echo "Not propper cell file format"
     exit 1
   fi

   #zcat $longfn >> $OUTPUT_DIR/$merged_file
   cat $longfn >> $OUTPUT_DIR/$merged_file.gz
   echo "cat $fn >> $merged_file.gz" | tee -a $LOG_FILE
   #zcat $longfn | cut -c -75 | gzip -c > $OUTPUT_DIR/$fn
done
