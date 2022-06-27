#!/bin/bash
name="MY_COHORT"
gtf="/mgh/references/gencode.v34.annotation.gtf"
jar="/mgh/gCNV/gCNV_helper.jar"
memb="-Xmx30G"
java ${memb} -jar ${jar} convertVCFsToBEDFormat -d ./ -o ${name}_vcftobed.bed -p genotyped-segments- -s .vcf
java ${memb} -jar ${jar} bedcluster -i ${name}_vcftobed.bed  -o ${name}_bedclustered.bed
java ${memb} -jar ${jar} defragment -i ${name}_bedclustered.bed -o ${name}_defragmented.bed filtered_intervals.txt 
java ${memb} -jar ${jar} calculateFrequency -i ${name}_defragmented.bed -c variant_name -o ${name}_frequencied.bed
java ${memb} -jar ${jar} updateNP -i ${name}_frequencied.bed -o ${name}_updateNP.bed -m filtered_intervals.txt
java ${memb} -jar ${jar} filter -i ${name}_updateNP.bed -o ${name}_filtered.bed
java ${memb} -jar ${jar} annotateWithGenes -a any -g ${gtf} -i ${name}_filtered.bed -o ${name}_intermediate1.bed -c svtype
java ${memb} -jar ${jar} countExons -i ${name}_intermediate1.bed -c genes_any_overlap -g ${gtf} -o ${name}_intermediate2.bed
java ${memb} -jar ${jar} annotateWithGenes -a strict -g ${gtf} -i ${name}_intermediate2.bed -o ${name}_intermediate3.bed -c svtype
java ${memb} -jar ${jar} countExons -i ${name}_intermediate3.bed -c genes_strict_overlap -g ${gtf} -o ${name}_annotated.bed
