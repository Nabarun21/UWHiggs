#!/bin/bash

# A stupid script for discovering the input ntuple files.
# Usage: ./get_inputs.sh jobid $hdfs
# Generates inputs/[jobid]/[sample1].txt
#           inputs/[jobid]/[sample2].txt
#           etc

JOBID=$1 
SOURCE=$2

echo "Finding input files for job: $JOBID in $SOURCE"

OUTPUTDIR=inputs/$JOBID

mkdir -p $OUTPUTDIR

for sampledir in `ls -d $SOURCE/$JOBID/*`
do
  sample=`basename $sampledir`
  echo -n "Getting data files for $sample - got "
  ls $sampledir/*/*root | sed "s|^|file:|" | root_file_check.py > $OUTPUTDIR/${sample}.tmp.txt
  echo -n `cat $OUTPUTDIR/${sample}.tmp.txt | wc -l `
  echo -n " files, "
  # Only update the file if it's changed
  cmp -s $OUTPUTDIR/${sample}.tmp.txt $OUTPUTDIR/${sample}.txt 
  if [ $? -ne 0 ]; then
    echo "new files found."
    mv $OUTPUTDIR/${sample}.tmp.txt $OUTPUTDIR/${sample}.txt 
  else
    echo "no new files found."
    rm $OUTPUTDIR/${sample}.tmp.txt 
  fi

done
