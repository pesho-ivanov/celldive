#!/bin/bash
# Runs tracer on all the reads

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
CONFIG_FILE=/home/pesho/libs/tracer-0.4.0/tracer-human.conf
CORES=20

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir $OUTPUT_DIR

echo "### run-human.sh ###"
echo "Processing all .fastq(.gz) files in ${INPUT_DIR}"
echo "Output in ${OUTPUT_DIR}"

# tracer python venv
#source ~/libs/tracer-0.4.0/bin/activate

tracer=~/libs/tracer-0.4.0/tracer

COUNTER=0
# Run tracer assemble for every cell
for longfn in $INPUT_DIR/*.fastq.gz
do
   fn=${longfn##*/}
   base=${fn%.fastq.gz}
   if [ ${base:(-1)} == "1" ]; then
     echo "`date`"
     echo $longfn ${longfn/_1.fastq.gz/_2.fastq.gz}
     PATH=$PATH:/home/pesho/libs/bowtie-1.1.2 IGDATA=//home/pesho/libs/igblast-1.4.0/ $tracer assemble -c $CONFIG_FILE -p 1 -s Hsap $longfn ${longfn/1.fastq.gz/2.fastq.gz} $base $OUTPUT_DIR 2>&1 | tee $LOG_FILE &
     let COUNTER=COUNTER+1
     if (( $COUNTER % $CORES == 0 )); then wait; fi
   fi
done

# waiting for all parallel cells
wait

# Summarise
$tracer summarise -c $CONFIG_FILE $OUTPUT_DIR

