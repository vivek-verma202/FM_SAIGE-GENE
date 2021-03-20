#!/bin/bash
# eg: /home/xao1/project/Rare_variant_pain/script/run_saige_step1.sh binary BackPain
# eg: /home/xao1/project/Rare_variant_pain/script/run_saige_step1.sh quantitative multisite
 
docker run --user $(id -u):$(getent group ukb | cut -d ':' -f 3) \
     	-v /scratch/xao1/project/Rare_variant_pain/saige_50k:/mnt/plink \
     	-v /scratch/xao1/project/Rare_variant_pain:/mnt/ \
     	-w /scratch/xao1/project/Rare_variant_pain/saige_50k \
		6eef0eeb60e2 nice -n 10 Rscript /usr/local/bin/step1_fitNULLGLMM.R     \
		--plinkFile=/mnt/plink/01.efe_caucasian \
        --phenoFile=/mnt/pheno/covars_and_pheno_${2}_chrnc.txt \
        --phenoCol=PHENO \
        --covarColList=SEX,AGE,AGE_SQ,UKBL,PCA1,PCA2,PCA3,PCA4,PCA5,PCA6,PCA7,PCA8,PCA9,PCA10,PCA11,PCA12,PCA13,PCA14,PCA15,PCA16,PCA17,PCA18,PCA19,PCA20,PCA21,PCA22,PCA23,PCA24,PCA25,PCA26,PCA27,PCA28,PCA29,PCA30,PCA31,PCA32,PCA33,PCA34,PCA35,PCA36,PCA37,PCA38,PCA39,PCA40,C_10003,C_11001,C_11002,C_11003,C_11004,C_11005,C_11006,C_11007,C_11008,C_11009,C_11011,C_11012,C_11013,C_11014,C_11016,C_11017,C_11018,C_11020,C_11021,C_11022,C_11023 \
        --sampleIDColinphenoFile=IID \
        --traitType=${1}       \
        --invNormalize=TRUE     \
        --outputPrefix=/mnt/plink/00_output/00_${2}_chrnc/out_${2}_chrnc \
		--outputPrefix_varRatio=/mnt/plink/00_output/00_${2}_chrnc/out_${2}_chrnc_cate \
		--sparseGRMFile=/mnt/plink/00_output/efe_caucasian_sparseGRM_numRandomMarkerforSparseKin_2000_relatednessCutoff_0.125_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx    \
        --sparseGRMSampleIDFile=/mnt/plink/00_output/efe_caucasian_sparseGRM_numRandomMarkerforSparseKin_2000_relatednessCutoff_0.125_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
        --nThreads=27 \
        --LOCO=FALSE \
		--skipModelFitting=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE	
