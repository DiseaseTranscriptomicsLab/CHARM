#!/bin/bash
rmats.py --b2 samples_unstressed.txt \
      --b1 samples_stressed.txt \
			--gtf ~/Genomes/Human/gencode.v42.primary_assembly.basic.annotation.gtf \
			--bi ~/Genomes/Human/STARindex/ \
      -t single \
			--readLength 76 \
			--nthread 4 \
			--od final_output_rMATS \
			--tmp intermediate_output_rMATS > rMATS_GSE171009.txt 