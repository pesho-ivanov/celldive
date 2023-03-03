#!/bin/bash

if [ $# -ne 2 ]; then
  echo 'Please give two arguments:'
  echo '  -- input dir'
  echo '  -- output dir'
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

if [ -d "$OUTPUT_DIR"  ]; then
  echo "The output dir $OUTPUT_DIR already exists. Please delete or change"
  exit 1
fi

mkdir -p $OUTPUT_DIR

echo "### trimmomatic.sh ###"
echo "Processing all .fastq(.gz) files in ${INPUT_DIR}"
echo "Output in ${OUTPUT_DIR}"

COUNTER=0
# Run tracer assemble for every cell
for longfn in $INPUT_DIR/*.fastq*
do
    fn=${longfn##*/}
    firstbase=${fn%.fastq.gz}
    if [ ${firstbase:(-2)} == "_1" ];
    then
        base=${firstbase::(-2)}
        second=${longfn/_1.fastq.gz/_2.fastq.gz}
        echo $longfn $second
        mkdir $OUTPUT_DIR/PE
        mkdir $OUTPUT_DIR/SE
        java -jar $libs/trimmomatic-0.36/trimmomatic-0.36.jar \
            PE -threads 20 $longfn $second \
            $OUTPUT_DIR/PE/${base}_1.fastq.gz $OUTPUT_DIR/SE/${base}_1.fastq.gz \
            $OUTPUT_DIR/PE/${base}_2.fastq.gz $OUTPUT_DIR/SE/${base}_2.fastq.gz \
            ILLUMINACLIP:$libs/trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 \
            MINLEN:100
   fi
done

# PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
