#!/bin/bash

# Isaac Wong
# iwong1@mgh.harvard.edu
# Version 1.31 February 21, 2020
# https://github.com/theisaacwong/

date

if [ "$1" == "-h" ]; then
  echo "arg 1: fastq file 1"
  echo "arg 2: fastq file 2"
  echo "arg 3: name of sample"
  echo "arg 4: path to output directory"
  echo "arg 5: path to a temp directory for storing intermediate files"
  echo "arg 6: amount of memory to give to JVM"
  exit 0
fi


#######################################
############ ENV VARIABLES ############
#######################################

FASTQ_PATH_R1=$1
FASTQ_PATH_R2=$2
NAME=$3
TMP_DIR=$4
TMP=$5
MEM=$6
UMI1="3M2S+T" 
UMI2="3M2S+T"
GOTO_STEP=${7:-0}

cd $TMP_DIR


NT=$(nproc)
PARIFQ="/data/talkowski/iwong/software/pairfq_lite/pairfq_lite"
REF="/data/talkowski/iwong/files/broad_references/Homo_sapiens_assembly38.fasta"
PICARD_PATH="/data/talkowski/iwong/src/jars/picard.jar"
FGBIO_PATH="/data/talkowski/iwong/software/fgbio/fgbio-1.0.0.jar"
GATK_PATH="/data/talkowski/iwong/src/jars/gatk.jar"
BAIT="/data/talkowski/iwong/files/intervals/Twist_Exome_Target_hg38.interval_list"
GNOMAD="/data/talkowski/iwong/files/references/gnomad.genomes.r3.0.sites.vcf.bgz"

LOG_DIR="${TMP_DIR}log_files/"
mkdir "${LOG_DIR}"

BSUB_ARGS=(-q big -n 1 -sla miket_sc -R 'rusage[swp=16000]' -R 'rusage[mem=64000]' -J ${NAME}_metrics -o ${LOG_DIR}metrics)
JAVA_ARGS=("-Dpicard.useLegacyParser=false" "-Djava.io.tmpdir=${TMP}" -Xmx${MEM} -jar)
FASTQC="/data/talkowski/iwong/software/FastQC/fastqc"

#######################################
########### OUPUT VARIABLES ###########
#######################################
READ_GROUP=${NAME}
STEP_02_OUTPUT=${TMP_DIR}${NAME}_step02_unaligned.bam
STEP_03_OUTPUT=${TMP_DIR}${NAME}_step03_umi_extracted_unaligned.bam
STEP_04_OUTPUT=${TMP_DIR}${NAME}_step04_queryname_sorted.ubam
STEP_05_OUTPUT=${TMP_DIR}${NAME}_step05_markAdapt.ubam
STEP_06_OUTPUT_01=${TMP_DIR}${NAME}_step06_R1.fastq
STEP_06_OUTPUT_02=${TMP_DIR}${NAME}_step06_R2.fastq
STEP_07_OUTPUT=${TMP_DIR}${NAME}_step07_mapped.sam
STEP_08_OUTPUT=${TMP_DIR}${NAME}_step08_merged.bam
STEP_09_OUTPUT=${TMP_DIR}${NAME}_step09_groupByUMI.bam
STEP_16_OUTPUT=${TMP_DIR}${NAME}_CollectDuplexSeqMetrics.txt
STEP_10_OUTPUT=${TMP_DIR}${NAME}_step10_callDuplexConsensus.bam
STEP_11_OUTPUT_01=${TMP_DIR}${NAME}_step11_R1.fastq
STEP_11_OUTPUT_02=${TMP_DIR}${NAME}_step11_R2.fastq
STEP_12_OUTPUT=${TMP_DIR}${NAME}_step12_mapped.sam
STEP_13_OUTPUT=${TMP_DIR}${NAME}_step13_merged.bam
STEP_14_OUTPUT=${TMP_DIR}${NAME}_step14_filterConsensus.bam
STEP_15_OUTPUT=${TMP_DIR}${NAME}_step15_clipBam.bam
STEP_16_OUTPUT=${TMP_DIR}${NAME}_step16_HaplotypeCaller.output.g.vcf.gz

#######################################
########### METRIC VARIABLES ##########
#######################################
STEP_02_METRICS_01=${TMP_DIR}${NAME}_step02_est_lib_complex_metrics.txt
STEP_03_METRICS_01=${TMP_DIR}${NAME}_step03_est_lib_complex_metrics.txt
STEP_03_METRICS_01=${TMP_DIR}${NAME}_step03_est_lib_complex_metrics.txt
STEP_08_METRICS_01=${TMP_DIR}${NAME}_step08_marked_duplicates.bam
STEP_08_METRICS_02=${TMP_DIR}${NAME}_step08_marked_dup_metrics.txt
STEP_08_METRICS_03=${TMP_DIR}${NAME}_step08_est_lib_complex_metrics.txt
STEP_08_METRICS_04=${TMP_DIR}${NAME}_step08_insert_size_metrics.txt
STEP_08_METRICS_05=${TMP_DIR}${NAME}_step08_insert_size_histogram.pdf
STEP_08_METRICS_06=${TMP_DIR}${NAME}_step08_hs_metrics.txt
STEP_10_METRICS_01=${TMP_DIR}${NAME}_step10_est_lib_complex_metrics.txt
STEP_13_METRICS_01=${TMP_DIR}${NAME}_step13_marked_duplicates.bam
STEP_13_METRICS_02=${TMP_DIR}${NAME}_step13_marked_dup_metrics.txt
STEP_13_METRICS_03=${TMP_DIR}${NAME}_step13_est_lib_complex_metrics.txt
STEP_13_METRICS_04=${TMP_DIR}${NAME}_step13_insert_size_metrics.txt
STEP_13_METRICS_05=${TMP_DIR}${NAME}_step13_insert_size_histogram.pdf
STEP_13_METRICS_06=${TMP_DIR}${NAME}_step13_hs_metrics.txt
STEP_14_METRICS_01=${TMP_DIR}${NAME}_step14_hs_metrics.txt
STEP_15_METRICS_01=${TMP_DIR}${NAME}_step15_est_lib_complex_metrics.txt
STEP_15_METRICS_02=${TMP_DIR}${NAME}_step15_hs_metrics.txt
STEP_15_METRICS_03=${TMP_DIR}${NAME}_step15_CollectInsertSizeMetrics.txt
STEP_15_METRICS_04=${TMP_DIR}${NAME}_step15_CollectInsertSizeMetrics.pdf



echo "Step 2"; date	
if [ $GOTO_STEP -le 2 ]; then
# Step 2, convert R1 and R2 to ubam
java "${JAVA_ARGS[@]}" ${PICARD_PATH} FastqToSam -O ${STEP_02_OUTPUT} -F1 ${FASTQ_PATH_R1} -F2 ${FASTQ_PATH_R2} -SM ${NAME} -LB Library1 -PU Unit1 -PL Illumina -RG ${READ_GROUP}
fi

echo "Step 3"; date
if [ $GOTO_STEP -le 3 ]; then
# Step 3 FGBIO ExtractUmisFromBam
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${FGBIO_PATH}" ExtractUmisFromBam --input=${STEP_02_OUTPUT} --output=${STEP_03_OUTPUT} --read-structure=${UMI1} ${UMI2} --molecular-index-tags=ZA ZB --single-tag=RX
fi

echo "Step 4"; date
if [ $GOTO_STEP -le 4 ]; then
# Step 4 sort ubam by query name
java "${JAVA_ARGS[@]}" "${PICARD_PATH}" SortSam -I ${STEP_03_OUTPUT} -O ${STEP_04_OUTPUT} -SORT_ORDER queryname
fi

echo "Step 5"; date
if [ $GOTO_STEP -le 5 ]; then
# Step 5 Mark Illumina Adapters
java "${JAVA_ARGS[@]}" "${PICARD_PATH}" MarkIlluminaAdapters -I ${STEP_04_OUTPUT} -O ${STEP_05_OUTPUT} -M ${STEP_05_OUTPUT}_MarkAdaptMetrics.txt
fi

echo "Step 6"; date
if [ $GOTO_STEP -le 6 ]; then
# Step 6, Sam to fastq
java "${JAVA_ARGS[@]}" "${PICARD_PATH}" SamToFastq -I ${STEP_05_OUTPUT} -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION X -CLIPPING_MIN_LENGTH 36 -INCLUDE_NON_PF_READS true -F ${STEP_06_OUTPUT_01} -F2 ${STEP_06_OUTPUT_02}
fi

echo "Step 7"; date
if [ $GOTO_STEP -le 7 ]; then
# Step 7, map fastq 
bwa mem -C -t ${NT} ${REF} ${STEP_06_OUTPUT_01} ${STEP_06_OUTPUT_02} > ${STEP_07_OUTPUT}
fi

echo "Step 8"; date		
if [ $GOTO_STEP -le 8 ]; then
CONTINUE_ON=true	
# Step 8, merge bam alignment
java "${JAVA_ARGS[@]}" "${PICARD_PATH}" MergeBamAlignment -UNMAPPED ${STEP_03_OUTPUT} -ALIGNED ${STEP_07_OUTPUT} -O ${STEP_08_OUTPUT} -R ${REF} -CLIP_ADAPTERS false -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true -EXPECTED_ORIENTATIONS FR -MAX_GAPS -1 -SO coordinate -ALIGNER_PROPER_PAIR_FLAGS false
fi


echo "Step 9"; date				# if you want a fastq files, after this step is the time to make them
if [ $GOTO_STEP -le 9 ]; then
# Step 9, Group Reads by UMI
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${FGBIO_PATH}" GroupReadsByUmi --strategy=paired --input=${STEP_08_OUTPUT} --output=${STEP_09_OUTPUT} --raw-tag=RX --assign-tag=MI --min-map-q=10 --edits=1
fi

echo "Step 10"; date
if [ $GOTO_STEP -le 10 ]; then
# Step 10, Call duplex consensus reads
#java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar /data/talkowski/iwong/scripts/samToFastqKeepUMI_jre1.8.jar "${TMP_DIR}${NAME}_groupedUMI_step07.bam" "${TMP_DIR}${NAME}_groupedUMI_step07.5.fastq"
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${FGBIO_PATH}" CallDuplexConsensusReads --input=${STEP_09_OUTPUT} --output=${STEP_10_OUTPUT} --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=10 --threads=${NT} --min-reads=0
fi

echo "Step 11"; date
if [ $GOTO_STEP -le 11 ]; then
# Step 11, remap duplex consensus reads
java "${JAVA_ARGS[@]}" "${PICARD_PATH}" SamToFastq -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true -I ${STEP_10_OUTPUT} -F ${STEP_11_OUTPUT_01} -F2 ${STEP_11_OUTPUT_02} -INCLUDE_NON_PF_READS true
fi

echo "Step 12"; date
if [ $GOTO_STEP -le 12 ]; then
# Step 12, map fastq 
bwa mem -C -t ${NT} ${REF} ${STEP_11_OUTPUT_01} ${STEP_11_OUTPUT_02} > ${STEP_12_OUTPUT}
fi

echo "Step 13"; date
if [ $GOTO_STEP -le 13 ]; then
# Step 13, merge bam alignment
java "${JAVA_ARGS[@]}" ${PICARD_PATH} MergeBamAlignment -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true -UNMAPPED ${STEP_10_OUTPUT} -ALIGNED ${STEP_12_OUTPUT} -OUTPUT ${STEP_13_OUTPUT} -REFERENCE_SEQUENCE ${REF} -CLIP_ADAPTERS false -ORIENTATIONS FR -MAX_GAPS -1 -SORT_ORDER coordinate -ALIGNER_PROPER_PAIR_FLAGS false -ATTRIBUTES_TO_RETAIN X0 -ATTRIBUTES_TO_RETAIN ZS -ATTRIBUTES_TO_RETAIN ZI -ATTRIBUTES_TO_RETAIN ZM -ATTRIBUTES_TO_RETAIN ZC -ATTRIBUTES_TO_RETAIN ZN -ATTRIBUTES_TO_REVERSE ad -ATTRIBUTES_TO_REVERSE bd -ATTRIBUTES_TO_REVERSE cd -ATTRIBUTES_TO_REVERSE ae -ATTRIBUTES_TO_REVERSE be -ATTRIBUTES_TO_REVERSE ce
fi


echo "Step 14"; date
if [ $GOTO_STEP -le 14 ]; then
# Step 14, FilterConsensusReads
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${FGBIO_PATH}" FilterConsensusReads -i ${STEP_13_OUTPUT} -o ${STEP_14_OUTPUT} --ref ${REF} --min-reads 0 0 0 --max-read-error-rate 0.35 --max-base-error-rate 0.3 --min-base-quality 40 --max-no-call-fraction 0.25
fi

echo "Step 15"; date
if [ $GOTO_STEP -le 15 ]; then
# Step 15, clip bam files # will need an alt file with no clipping of overlapping-reads for variant calling
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${FGBIO_PATH}" ClipBam --input=${STEP_14_OUTPUT} --output=${STEP_15_OUTPUT} --ref=${REF} --clipping-mode=Hard --clip-overlapping-reads=true
bsub "${BSUB_ARGS[@]}" java "${JAVA_ARGS[@]}" "${PICARD_PATH}" EstimateLibraryComplexity -I ${STEP_15_OUTPUT} -O ${STEP_15_METRICS_01}
bsub "${BSUB_ARGS[@]}" java "${JAVA_ARGS[@]}" "${PICARD_PATH}" CollectHsMetrics -I ${STEP_15_OUTPUT} -O ${STEP_15_METRICS_02} -R ${REF} -BAIT_INTERVALS ${BAIT} -TARGET_INTERVALS ${BAIT}
bsub "${BSUB_ARGS[@]}" java "${JAVA_ARGS[@]}" "${PICARD_PATH}" CollectInsertSizeMetrics -I ${STEP_15_OUTPUT} -O ${STEP_15_METRICS_03} -H ${STEP_15_METRICS_04}
fi
date

echo "Step 16"; date
if [ $GOTO_STEP -le 16 ]; then
# Step 15, clip bam files # will need an alt file with no clipping of overlapping-reads for variant calling
java -Xmx${MEM} -Djava.io.tmpdir=${TMP} -jar "${GATK_PATH}" HaplotypeCaller -R ${REF} -I ${STEP_15_OUTPUT} -O ${STEP_16_OUTPUT} -ERC GVCF --dbsnp ${GNOMAD} -L ${BAIT}
fi
date



