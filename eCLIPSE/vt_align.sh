#!/bin/bash


if [ -z "$1" ]
  then echo "No target directory has been provided. Assuming current path!"
  path=./
else
  path="$1"
fi

samples=$(ls ${path}*.paired.fastq.gz | sed  's/_1.paired.fastq.gz$//' | sed  's/_2.paired.fastq.gz$//' | sort -u | xargs -n 1 basename)

for sample in ${samples[@]}
do
	echo ${sample}
done

if [ ! -d ${path}/alignment_VASTTOOLS ]
	then mkdir ${path}/alignment_VASTTOOLS
fi

# Standard Unique mode
for sample in ${samples[@]}
do
  ulimit -n 10000
  vast-tools align ${path}/${sample}_1.paired.fastq.gz ${path}/${sample}_2.paired.fastq.gz \
 	--sp hg38 \
 	--name ${sample} \
 	--expr \
 	-c 60 > VT_align_${sample}.txt 2>&1
done