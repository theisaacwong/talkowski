#!/bin/bash

module load R/3.5.1-foss-2018b

WD="/data/talkowski/iwong/SFARI/"

cd ${WD}"xhmm"

XHMM="/PHShome/iw068/statgen-xhmm-cc14e528d909"
DATA_RD_txt="./DATA.RD.txt"
DATA_filtered_centered_RD_txt="./DATA.filtered_centered.RD.txt"
DATA_filtered_centered_RD_txt_filtered_targets_txt="./DATA.filtered_centered.RD.txt.filtered_targets.txt"
DATA_filtered_centered_RD_txt_filtered_samples_txt="./DATA.filtered_centered.RD.txt.filtered_samples.txt"
DATA_RD_PCA="./DATA.RD_PCA"
DATA_PCA_normalized_txt="./DATA.PCA_normalized.txt"
DATA_PCA_normalized_filtered_sample_zscores_RD_txt="./DATA.PCA_normalized.filtered.sample_zscores.RD.txt"
DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_targets_txt="./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"
DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_samples_txt="./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"
DATA_same_filtered_RD_txt="./DATA.same_filtered.RD.txt"
params_txt="./params.txt"
DATA_xcnv="./DATA.xcnv"
DATA_aux_xcnv="./DATA.aux_xcnv"
DATA="./DATA"
DATA_vcf="./DATA.vcf"



# this step converts the counts file into the matrix input for xhmm
#Rscript /data/talkowski/iwong/SFARI/xhmm/xhmm_prep.R /data/talowski/iwong/SFARI/xhmm/counts/ $DATA_RD_txt
#sed -i '1s;^;GATK._mean_cvg;' $DATA_RD_txt


# Filters samples and targets and then mean-centers the targets:
$XHMM/xhmm --matrix -r $DATA_RD_txt --centerData --centerType target \
-o $DATA_filtered_centered_RD_txt \
--outputExcludedTargets $DATA_filtered_centered_RD_txt_filtered_targets_txt \
--outputExcludedSamples $DATA_filtered_centered_RD_txt_filtered_samples_txt \
--minTargetSize 10 --maxTargetSize 10000  \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150


# Runs PCA on mean-centered data:
$XHMM/xhmm --PCA -r $DATA_filtered_centered_RD_txt --PCAfiles $DATA_RD_PCA


# Normalizes mean-centered data using PCA information:
$XHMM/xhmm --normalize -r $DATA_filtered_centered_RD_txt --PCAfiles $DATA_RD_PCA \
--normalizeOutput $DATA_PCA_normalized_txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7


# Filters and z-score centers (by sample) the PCA-normalized data:
$XHMM/xhmm --matrix -r $DATA_PCA_normalized_txt --centerData --centerType sample --zScoreData \
-o $DATA_PCA_normalized_filtered_sample_zscores_RD_txt \
--outputExcludedTargets $DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_targets_txt \
--outputExcludedSamples $DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_samples_txt \
--maxSdTargetRD 30


# Filters original read-depth data to be the same as filtered, normalized data:
$XHMM/xhmm --matrix -r $DATA_RD_txt \
--excludeTargets $DATA_filtered_centered_RD_txt_filtered_targets_txt \
--excludeTargets $DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_targets_txt \
--excludeSamples $DATA_filtered_centered_RD_txt_filtered_samples_txt \
--excludeSamples $DATA_PCA_normalized_filtered_sample_zscores_RD_txt_filtered_samples_txt \
-o $DATA_same_filtered_RD_txt


# Discovers CNVs in normalized data:
$XHMM/xhmm --discover -p $params_txt \
-r $DATA_PCA_normalized_filtered_sample_zscores_RD_txt -R $DATA_same_filtered_RD_txt \
-c $DATA_xcnv -a $DATA_aux_xcnv -s $DATA


# NOTE: A description of the .xcnv file format can be found below.


# Genotypes discovered CNVs in all samples:
$XHMM/xhmm --genotype -p $params_txt \
-r $DATA_PCA_normalized_filtered_sample_zscores_RD_txt -R $DATA_same_filtered_RD_txt \
-g $DATA_xcnv -F /data/talkowski/iwong/hg19/Homo_sapiens_assembly19.fasta \
-v $DATA_vcf





