#!/bin/bash

set -eou pipefail

name="MGH"
gtf="/mgh/references/gencode.v34.annotation.gtf"
jar="/mgh/scripts/gCNV_helper.jar"
memb="-Xmx30G"
java ${memb} -jar ${jar} downloadSegmentsVCFs -m sample_set_entity.tsv -d ./ -c genotyped_segments_vcfs
java ${memb} -jar ${jar} downloadFilteredIntervals -m sample_set_entity.tsv -d ./ -c filtered_intervals -o filtered_intervals.txt
java ${memb} -jar ${jar} convertVCFsToBEDFormat -d ./ -o ${name}_vcftobed.bed -p genotyped-segments- -s .vcf
java ${memb} -jar ${jar} defragment -i ${name}_vcftobed.bed -o ${name}_defragmented.bed -m filtered_intervals.txt
java ${memb} -jar ${jar} bedcluster -i ${name}_defragmented.bed  -o ${name}_bedclustered.bed
java ${memb} -jar ${jar} calculateFrequency -i ${name}_bedclustered.bed -c variant_name -o ${name}_frequencied.bed
java ${memb} -jar ${jar} updateNP -i ${name}_frequencied.bed -o ${name}_updateNP.bed -m filtered_intervals.txt
java ${memb} -jar ${jar} filter -i ${name}_updateNP.bed -o ${name}_filtered.bed
java ${memb} -jar ${jar} annotateWithGenes -a any -g ${gtf} -i ${name}_filtered.bed -o ${name}_intermediate1.bed -c svtype
java ${memb} -jar ${jar} countExons -i ${name}_intermediate1.bed -c genes_any_overlap -g ${gtf} -o ${name}_intermediate2.bed
java ${memb} -jar ${jar} annotateWithGenes -a strict -g ${gtf} -i ${name}_intermediate2.bed -o ${name}_intermediate3.bed -c svtype
java ${memb} -jar ${jar} countExons -i ${name}_intermediate3.bed -c genes_strict_overlap -g ${gtf} -o ${name}_annotated.bed

rm -rf ${name}_intermediate1.bed
rm -rf ${name}_intermediate2.bed
rm -rf ${name}_intermediate3.bed

rm -rf ${name}_frequencied.bed
rm -rf ${name}_updateNP.bed
rm -rf ${name}_filtered.bed


