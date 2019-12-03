#!/bin/bash

date

if [ "$1" == "-h" ]; then
  echo "arg 1: fastq file 1"
  echo "arg 2: fastq file 2"
  echo "arg 3: name of sample"
  echo "arg 4: path to picard.jar"
  echo "arg 5: path to fgbio.jar"
  echo "arg 6: path to output directory"
  echo "arg 7: path to a temp directory for storing intermediate files"
  echo "arg 8: amount of memory to give to JVM"
  echo "arg 9: UMI flag"
  echo "arg 10: 2nd UMI flag (optional)"
  exit 0
fi


PICARD_PATH=$4
FGBIO_PATH=$5

#java -Xmx$MEM -jar $PICARD_PATH CreateSequenceDictionary R="${TMP_DIR}hg38.fa" O="${TMP_DIR}hg38.dict"

FASTQ_PATH_R1=$1
FASTQ_PATH_R2=$2

TMP_DIR=$6

TMP=$7

NAME=$3

MEM=$8

NT=$(nproc)

echo -e "\e[104mstep one, FastqToSam\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH FastqToSam O="${TMP_DIR}${NAME}_unaligned_step01.bam" F1=$FASTQ_PATH_R1 F2=$FASTQ_PATH_R2 SM=$NAME LB=Library1 PU=Unit1 PL=Illumina


echo -e "\e[104mstep two, ExtractUmisFromBam\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $FGBIO_PATH ExtractUmisFromBam --input="${TMP_DIR}${NAME}_unaligned_step01.bam" --output="${TMP_DIR}${NAME}_unaligned_umi_step02.bam" --read-structure=$9 ${10} --molecular-index-tags=ZA ZB --single-tag=RX


echo -e "\e[104mstep three, MarkIlluminaAdapters\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH MarkIlluminaAdapters I="${TMP_DIR}${NAME}_unaligned_umi_step02.bam" O="${TMP_DIR}${NAME}_unaligned_umi_markedAdpt_step03.bam" M="${TMP_DIR}${NAME}_adapter_Metrics_step03.txt"


echo -e "\e[104mstep four, SamToFastq\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH SamToFastq I="${TMP_DIR}${NAME}_unaligned_umi_markedAdpt_step03.bam" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X CLIPPING_MIN_LENGTH=36 INCLUDE_NON_PF_READS=true INTERLEAVE=true F="${TMP_DIR}${NAME}_unaligned_umi_markedAdpt_step04.fastq"


echo -e "\e[104mstep five, bwa mem 1\e[0m"
	bwa mem -p -t $NT "${TMP_DIR}hg38.fa" "${TMP_DIR}${NAME}_unaligned_umi_markedAdpt_step04.fastq" > "${TMP_DIR}${NAME}_aligned_umi_markedAdpt_bwa_mem_step05.sam"


echo -e "\e[104mstep six, MergeBamAlignment\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH MergeBamAlignment UNMAPPED="${TMP_DIR}${NAME}_unaligned_umi_step02.bam" ALIGNED="${TMP_DIR}${NAME}_aligned_umi_markedAdpt_bwa_mem_step05.sam" O="${TMP_DIR}${NAME}_mapped_step06.bam" R="${TMP_DIR}hg38.fa" CLIP_ADAPTERS=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true EXPECTED_ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false


echo -e "\e[104mstep seven, GroupReadsByUmi\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $FGBIO_PATH GroupReadsByUmi --strategy=paired --input="${TMP_DIR}${NAME}_mapped_step06.bam" --output="${TMP_DIR}${NAME}_groupedUMI_step07.bam" --raw-tag=RX --assign-tag=MI --min-map-q=10 --edits=1

echo -e "\e[104mstep seven 1/2 \e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar /data/talkowski/iwong/scripts/samToFastqKeepUMI_jre1.8.jar "${TMP_DIR}${NAME}_groupedUMI_step07.bam" "${TMP_DIR}${NAME}_groupedUMI_step07.5.fastq"
	
#echo -e "\e[104mstep eight, CallDuplexConsensusReads\e[0m" 
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $FGBIO_PATH CallDuplexConsensusReads --input="${TMP_DIR}${NAME}_groupedUMI_step07.bam" --output="${TMP_DIR}${NAME}_ds_consensus_unaligned_step08.bam" --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=10 --threads=${NT} --min-reads=0


#echo -e "\e[104mstep nine,SamToFastq\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH SamToFastq VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT="${TMP_DIR}${NAME}_ds_consensus_unaligned_step08.bam" F="${TMP_DIR}${NAME}_ds_consensus_unaligned_step09.fastq" INTERLEAVE=true INCLUDE_NON_PF_READS=true


#echo -e "\e[104mstep ten, bwa mem 2\e[0m"
	bwa mem -t $NT -p "${TMP_DIR}hg38.fa" "${TMP_DIR}${NAME}_ds_consensus_unaligned_step09.fastq" > "${TMP_DIR}${NAME}_ds_consensus_aligned_bwa_mem_step10.fastq"


#echo -e "\e[104mstep eleven, MergeBamAlignment\e[0m"
	java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar $PICARD_PATH MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true UNMAPPED="${TMP_DIR}${NAME}_ds_consensus_unaligned_step08.bam" ALIGNED="${TMP_DIR}${NAME}_ds_consensus_aligned_bwa_mem_step10.fastq" OUTPUT="${TMP_DIR}${NAME}_ds_consensus_aligned_step11.bam" REFERENCE_SEQUENCE="${TMP_DIR}hg38.fa" CLIP_ADAPTERS=false ORIENTATIONS=FR MAX_GAPS=-1 SORT_ORDER=coordinate ALIGNER_PROPER_PAIR_FLAGS=false ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN ATTRIBUTES_TO_REVERSE=ad ATTRIBUTES_TO_REVERSE=bd ATTRIBUTES_TO_REVERSE=cd ATTRIBUTES_TO_REVERSE=ae ATTRIBUTES_TO_REVERSE=be ATTRIBUTES_TO_REVERSE=ce


#echo -e "\e[104mstep twelve, FilterConsensusReads\e[0m"
	java -Xmx$MEM -jar $FGBIO_PATH FilterConsensusReads -i "${TMP_DIR}${NAME}_ds_consensus_aligned_step11.bam" -o "${TMP_DIR}${NAME}_ds_consensus_filtered_step12.bam" --ref "${TMP_DIR}hg38.fa" --min-reads 10 5 3 --max-read-error-rate 0.05 --max-base-error-rate 0.1 --min-base-quality 50 --max-no-call-fraction 0.05


#echo -e "\e[104mstep thirteen, ClipBam\e[0m"
	java -Xmx$MEM -jar $FGBIO_PATH ClipBam --input="${TMP_DIR}${NAME}_ds_consensus_filtered_step12.bam" --output="${TMP_DIR}${NAME}_ds_consensus_filtered_clipped_step13.bam" --ref="${TMP_DIR}hg38.fa" --clipping-mode=Hard --clip-overlapping-reads=true


#echo -e "\e[104mstep fourteen, CollectDuplexSeqMetrics\e[0m"
	java -Xmx$MEM -jar $FGBIO_PATH CollectDuplexSeqMetrics --input="${TMP_DIR}${NAME}_groupedUMI_step07.bam" --output="${TMP_DIR}${NAME}_duplexSeqMetrics" --description=$NAME

date


