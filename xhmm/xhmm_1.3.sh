#!/bin/bash
source activate iwong1

dir=/data/talkowski/iwong/SFARI/xhmm
xhmm=/PHShome/iw068/statgen-xhmm-cc14e528d909/xhmm
batch=OCT4

minMeanSampleRD=${1:-25}
maxMeanSampleRD=${2:-200}
maxSdSampleRD=${3:-150}
maxSdTargetRD=${4:-30}
PVE_mean_factor=${5:-0.7}
dir=${$6:-/data/talkowski/iwong/SFARI/xhmm}

echo minMeanSampleRD=${minMeanSampleRD}
echo maxMeanSampleRD=${maxMeanSampleRD}
echo maxSdSampleRD=${maxSdSampleRD}
echo maxSdTargetRD=${maxSdTargetRD}
echo PVE_mean_factor=${PVE_mean_factor}


################################################################################################################################################
## Filters samples and targets and then mean-centers the targets (filters the same targets as in train):
$xhmm --matrix -r ${dir}/${batch}.RD.txt --centerData --centerType target \
    -o ${dir}/${batch}.DATA.filtered_centered.RD.txt \
	--excludeTargets ${dir}/train.DATA.filtered_centered.RD.txt.filtered_targets.txt \
    --outputExcludedSamples ${dir}/${batch}.DATA.filtered_centered.RD.txt.filtered_samples.txt \
    --minMeanSampleRD ${minMeanSampleRD} --maxMeanSampleRD ${maxMeanSampleRD} \
    --maxSdSampleRD ${maxSdSampleRD}

################################################################################################################################################
# Runs PCA on mean-centered data:
$xhmm --PCA -r ${dir}/${batch}.DATA.filtered_centered.RD.txt --PCAfiles ${dir}/${batch}.RD.PCA.txt

################################################################################################################################################
## Normalizes mean-centered data using PCA information:
$xhmm --normalize -r ${dir}/${batch}.DATA.filtered_centered.RD.txt --PCAfiles ${dir}/${batch}.RD.PCA.txt \
    --normalizeOutput ${dir}/${batch}.DATA.PCA_normalized.txt \
    --PCnormalizeMethod PVE_mean --PVE_mean_factor ${PVE_mean_factor}

################################################################################################################################################
## Filters and z-score centers (by sample) the PCA-normalized data:
$xhmm --matrix -r ${dir}/${batch}.DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
    -o ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
    --outputExcludedTargets ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --outputExcludedSamples ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    --maxSdTargetRD ${maxSdTargetRD}

################################################################################################################################################
## Filters original read-depth data to be the same as filtered, normalized data:
$xhmm --matrix -r ${dir}/${batch}.RD.txt \
    --excludeTargets ${dir}/train.DATA.filtered_centered.RD.txt.filtered_targets.txt \
    --excludeTargets ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
    --excludeSamples ${dir}/${batch}.DATA.filtered_centered.RD.txt.filtered_samples.txt \
    --excludeSamples ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
    -o ${dir}/${batch}.DATA.same_filtered.RD.txt

################################################################################################################################################
## Discovers CNVs in normalized data:
$xhmm --discover -p ${dir}/params.txt \
    -r ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ${dir}/${batch}.DATA.same_filtered.RD.txt \
    -c ${dir}/${batch}.DATA.xcnv -a ${dir}/${batch}.DATA.aux_xcnv -s ${dir}/${batch}.DATA

################################################################################################################################################
## Genotypes discovered CNVs in all samples:
$xhmm --genotype -p ${dir}/params.txt \
    -r ${dir}/${batch}.DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
    -R ${dir}/${batch}.DATA.same_filtered.RD.txt \
    -g ${dir}/${batch}.DATA.xcnv \
    -F /data/talkowski/iwong/hg19/Homo_sapiens_assembly19.fasta \
    -v ${dir}/${batch}.DATA.vcf


