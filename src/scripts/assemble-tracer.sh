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
CONFIG_FILE=/home/pesho/libs/tracer-0.4.0/tracer-mouse.conf
CORES=15

LOG_FILE=$OUTPUT_DIR/tracer.log

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
   base=${fn%.fastq*}
   if [ ${base:(-2)} == "_1" ]; then
     echo "`date`"
     ./tracer assemble -c $CONFIG_FILE -p 1 -s Mmus $longfn ${longfn/_1.fastq/_2.fastq} $base $OUTPUT_DIR 2>&1 | tee $LOG_FILE &
     # source assemble-one.sh $longfn $base &
     let COUNTER=COUNTER+1
     if (( $COUNTER % $CORES == 0 )); then wait; fi
   fi 
done

# waiting for all parallel cells
wait

# Summarise
./tracer summarise -c $CONFIG_FILE $OUTPUT_DIR

