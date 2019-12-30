#!/bin/bash

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



FASTQ_PATH_R1=$1
FASTQ_PATH_R2=$2
TMP_DIR=$4
TMP=$5
NAME=$3
MEM=$6
UMI1="8M142T" 
UMI2="8M142T"

NT=$(nproc)
PARIFQ="/data/talkowski/iwong/software/pairfq_lite/pairfq_lite"
REF="/data/talkowski/iwong/files/hg38/hg38.fa"
PICARD_PATH="/data/talkowski/iwong/src/jars/picard.jar"
FGBIO_PATH="/data/talkowski/iwong/software/fgbio/fgbio-1.0.0.jar"


# Step 1, add /1 /2 suffix to R1/R2 PE reads
#echo "Step 1"
#STEP_01_OUTPUT_01=${TMP_DIR}R1_info.fq
#STEP_01_OUTPUT_02=${TMP_DIR}R2_info.fq
#sed 's, ,:,g' ${FASTQ_PATH_R1} > ${FASTQ_PATH_R1}_trm
#sed 's, ,:,g' ${FASTQ_PATH_R2} > ${FASTQ_PATH_R2}_trm
#${PARIFQ} addinfo -i ${FASTQ_PATH_R1}_trm -o ${STEP_01_OUTPUT_01} -p 1
#${PARIFQ} addinfo -i ${FASTQ_PATH_R2}_trm -o ${STEP_01_OUTPUT_02} -p 2

echo "Step 2"; date	
# Step 2, convert R1 and R2 to ubam
STEP_02_OUTPUT=${TMP_DIR}${NAME}_step02_unaligned.bam
STEP_02_METRICS_01=${TMP_DIR}${NAME}_step02_est_lib_complex_metrics.txt
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} FastqToSam O=${STEP_02_OUTPUT} F1=${FASTQ_PATH_R1} F2=${FASTQ_PATH_R2} SM=${NAME} LB=Library1 PU=Unit1 PL=Illumina
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_02_OUTPUT} O=${STEP_02_METRICS_01}

echo "Step 3"; date
# Step 3 FGBIO ExtractUmisFromBam
STEP_03_OUTPUT=${TMP_DIR}${NAME}_step03_umi_extracted_unaligned.bam
STEP_03_METRICS_01=${TMP_DIR}${NAME}_step03_est_lib_complex_metrics.txt
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${FGBIO_PATH} ExtractUmisFromBam --input=${STEP_02_OUTPUT} --output=${STEP_03_OUTPUT} --read-structure=${UMI1} ${UMI2} --molecular-index-tags=ZA ZB --single-tag=RX
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_03_OUTPUT} O=${STEP_03_METRICS_01}

echo "Step 4"; date
# Step 4 sort ubam by query name
STEP_04_OUTPUT=${TMP_DIR}${NAME}_step04_queryname_sorted.ubam
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} SortSam I=${STEP_03_OUTPUT} O=${STEP_04_OUTPUT} SORT_ORDER=queryname

echo "Step 5"; date
# Step 5 Mark Illumina Adapters
STEP_05_OUTPUT=${TMP_DIR}${NAME}_step05_markAdapt.ubam
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} MarkIlluminaAdapters I=${STEP_04_OUTPUT} O=${STEP_05_OUTPUT} M=${STEP_05_OUTPUT}_MarkAdaptMetrics.txt

echo "Step 6"; date
# Step 6, Sam to fastq
STEP_06_OUTPUT_01=${TMP_DIR}${NAME}_step06_R1.fastq
STEP_06_OUTPUT_02=${TMP_DIR}${NAME}_step06_R2.fastq
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} SamToFastq I=${STEP_05_OUTPUT} CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X CLIPPING_MIN_LENGTH=36 INCLUDE_NON_PF_READS=true F=${STEP_06_OUTPUT_01} F2=${STEP_06_OUTPUT_02}

echo "Step 7"; date
# Step 7, map fastq 
STEP_07_OUTPUT=${TMP_DIR}${NAME}_step07_mapped.sam
bwa mem -t ${NT} ${REF} ${STEP_06_OUTPUT_01} ${STEP_06_OUTPUT_02} > ${STEP_07_OUTPUT}

echo "Step 8"; date
# Step 8, merge bam alignment
STEP_08_OUTPUT=${TMP_DIR}${NAME}_step08_merged.bam
STEP_08_METRICS_01=${TMP_DIR}${NAME}_step08_marked_duplicates.bam
STEP_08_METRICS_02=${TMP_DIR}${NAME}_step08_marked_dup_metrics.txt
STEP_08_METRICS_03=${TMP_DIR}${NAME}_step08_est_lib_complex_metrics.txt
STEP_08_METRICS_04=${TMP_DIR}${NAME}_step08_insert_size_metrics.txt
STEP_08_METRICS_05=${TMP_DIR}${NAME}_step08_insert_size_histogram.pdf
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} MergeBamAlignment UNMAPPED=${STEP_03_OUTPUT} ALIGNED=${STEP_07_OUTPUT} O=${STEP_08_OUTPUT} R=${REF} CLIP_ADAPTERS=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true EXPECTED_ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} MarkDuplicates I=${STEP_08_OUTPUT} O=${STEP_08_METRICS_01} M=${STEP_08_METRICS_02}
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_08_OUTPUT} O=${STEP_08_METRICS_03}
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} CollectInsertSizeMetrics I=${STEP_08_OUTPUT} O=${STEP_08_METRICS_04} H=${STEP_08_METRICS_05}

echo "Step 9"; date
# Step 9, Group Reads by UMI
STEP_09_OUTPUT=${TMP_DIR}${NAME}_step09_groupByUMI.bam
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${FGBIO_PATH} GroupReadsByUmi --strategy=paired --input=${STEP_08_OUTPUT} --output=${STEP_09_OUTPUT} --raw-tag=RX --assign-tag=MI --min-map-q=10 --edits=1

echo "Step 10"; date
# Step 10, Call duplex consensus reads
#java -Djava.io.tmpdir=$TMP -Xmx$MEM -jar /data/talkowski/iwong/scripts/samToFastqKeepUMI_jre1.8.jar "${TMP_DIR}${NAME}_groupedUMI_step07.bam" "${TMP_DIR}${NAME}_groupedUMI_step07.5.fastq"
STEP_10_OUTPUT=${TMP_DIR}${NAME}_step10_callDuplexConsensus.bam
STEP_10_METRICS_01=${TMP_DIR}${NAME}_step10_est_lib_complex_metrics.txt
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${FGBIO_PATH} CallDuplexConsensusReads --input=${STEP_09_OUTPUT} --output=${STEP_10_OUTPUT} --error-rate-pre-umi=45 --error-rate-post-umi=30 --min-input-base-quality=10 --threads=${NT} --min-reads=0
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_10_OUTPUT} O=${STEP_10_METRICS_01}

echo "Step 11"; date
# Step 11, remap duplex consensus reads
STEP_11_OUTPUT_01=${TMP_DIR}${NAME}_step11_R1.fastq
STEP_11_OUTPUT_02=${TMP_DIR}${NAME}_step11_R2.fastq
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} SamToFastq VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true I=${STEP_10_OUTPUT} F=${STEP_11_OUTPUT_01} F2=${STEP_11_OUTPUT_02} INCLUDE_NON_PF_READS=true

echo "Step 12"; date
# Step 12, map fastq 
STEP_12_OUTPUT=${TMP_DIR}${NAME}_step12_mapped.sam
bwa mem -t ${NT} ${REF} ${STEP_11_OUTPUT_01} ${STEP_11_OUTPUT_02} > ${STEP_12_OUTPUT}

echo "Step 13"; date
# Step 13, merge bam alignment
STEP_13_OUTPUT=${TMP_DIR}${NAME}_step13_merged.bam
STEP_13_METRICS_01=${TMP_DIR}${NAME}_step13_marked_duplicates.bam
STEP_13_METRICS_02=${TMP_DIR}${NAME}_step13_marked_dup_metrics.txt
STEP_13_METRICS_03=${TMP_DIR}${NAME}_step13_est_lib_complex_metrics.txt
STEP_13_METRICS_04=${TMP_DIR}${NAME}_step13_insert_size_metrics.txt
STEP_13_METRICS_05=${TMP_DIR}${NAME}_step13_insert_size_histogram.pdf
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true UNMAPPED=${STEP_10_OUTPUT} ALIGNED=${STEP_12_OUTPUT} OUTPUT=${STEP_13_OUTPUT} REFERENCE_SEQUENCE=${REF} CLIP_ADAPTERS=false ORIENTATIONS=FR MAX_GAPS=-1 SORT_ORDER=coordinate ALIGNER_PROPER_PAIR_FLAGS=false ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN ATTRIBUTES_TO_REVERSE=ad ATTRIBUTES_TO_REVERSE=bd ATTRIBUTES_TO_REVERSE=cd ATTRIBUTES_TO_REVERSE=ae ATTRIBUTES_TO_REVERSE=be ATTRIBUTES_TO_REVERSE=ce
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} MarkDuplicates I=${STEP_13_OUTPUT} O=${STEP_13_METRICS_01} M=${STEP_13_METRICS_02}
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_13_OUTPUT} O=${STEP_13_METRICS_03}
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} CollectInsertSizeMetrics I=${STEP_13_OUTPUT} O=${STEP_13_METRICS_04} H=${STEP_13_METRICS_05}

echo "Step 14"; date
# Step 14, FilterConsensusReads
STEP_14_OUTPUT=${TMP_DIR}${NAME}_step14_filterConsensus.bam
java -Xmx${MEM} -jar ${FGBIO_PATH} FilterConsensusReads -i ${STEP_13_OUTPUT} -o ${STEP_14_OUTPUT} --ref ${REF} --min-reads 10 5 3 --max-read-error-rate 0.05 --max-base-error-rate 0.1 --min-base-quality 50 --max-no-call-fraction 0.05

echo "Step 15"; date
# Step 15, clip bam files 
STEP_15_OUTPUT=${TMP_DIR}${NAME}_step15_clipBam.bam
STEP_15_METRICS_01=${TMP_DIR}${NAME}_step15_est_lib_complex_metrics.txt
java -Xmx${MEM} -jar ${FGBIO_PATH} ClipBam --input=${STEP_14_OUTPUT} --output=${STEP_15_OUTPUT} --ref=${REF} --clipping-mode=Hard --clip-overlapping-reads=true
java -Djava.io.tmpdir=${TMP} -Xmx${MEM} -jar ${PICARD_PATH} EstimateLibraryComplexity I=${STEP_15_OUTPUT} O=${STEP_15_METRICS_01}

echo "Step 16"; date
# Step 16, CollectDuplexSeqMetrics
STEP_16_OUTPUT=${TMP_DIR}${NAME}_CollectDuplexSeqMetrics.txt
java -Xmx${MEM} -jar ${FGBIO_PATH} CollectDuplexSeqMetrics --input=${STEP_09_OUTPUT} --output=${STEP_16_OUTPUT} --description=${NAME}

date


