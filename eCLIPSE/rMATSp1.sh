#!/bin/bash

if [ -z "$1" ]
  then echo "No target directory has been provided. Assuming current path!"
  path=./
else
  path="$1"
fi

samples=$(ls ${path}*.bam | sed  's/.bam$//' | sort -u | xargs -n 1 basename)

for sample in ${samples[@]}
do
	echo ${sample}
done

for sample in ${samples[@]}
do samtools index -b ${path}/${sample}.bam
done