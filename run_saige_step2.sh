#!/bin/bash

# # test saige data
# docker run --user $(id -u):$(getent group ukb | cut -d ':' -f 3) \
#         -v /scratch/xao1/project/Rare_variant_pain/saige_50k:/mnt/plink \
#         -v /scratch/xao1/project/Rare_variant_pain:/mnt/ \
#         -v /scratch/xao1/tools/SAIGE/extdata:/mnt/test/ \
#         -w /scratch/xao1/project/Rare_variant_pain/saige_50k \
#         6eef0eeb60e2 nice -n 10 Rscript /usr/local/bin/step2_SPAtests.R \
#         --vcfFile=/mnt/test/input/genotype_10markers.vcf.gz \
#         --vcfFileIndex=/mnt/test/input/genotype_10markers.vcf.gz.csi \
#         --vcfField=GT \
#         --chrom=1 \
#         --minMAF=0 \
#         --minMAC=0.5 \
#         --maxMAFforGroupTest=0.01       \
#         --sampleFile=/mnt/test/input/samplelist.txt \
#         --GMMATmodelFile=/mnt/test/output/example_quantitative.rda \
#         --varianceRatioFile=/mnt/test/output/example_quantitative_cate.varianceRatio.txt \
#         --SAIGEOutputFile=/mnt/plink/test/example_quantitative.SAIGE.gene.txt \
#         --numLinesOutput=1 \
#         --groupFile=/mnt/test/input/groupFile_geneBasedtest_simulation_edited.txt    \
#         --sparseSigmaFile=/mnt/test/output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
#         --IsSingleVarinGroupTest=TRUE \
#         --LOCO=FALSE

# my run
docker run --user $(id -u):$(getent group ukb | cut -d ':' -f 3) \
        -v /scratch/xao1/project/Rare_variant_pain/saige_50k:/mnt/plink \
        -v /scratch/xao1/project/Rare_variant_pain:/mnt/ \
        -w /scratch/xao1/project/Rare_variant_pain/saige_50k \
        6eef0eeb60e2 nice -n 10 Rscript /usr/local/bin/step2_SPAtests.R \
        --vcfFile=/mnt/plink/01.efe_caucasian_${2}.vcf.gz \
        --vcfFileIndex=/mnt/plink/01.efe_caucasian_${2}.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=${2} \
        --minMAF=0 \
        --minMAC=3 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=/mnt/plink/01.efe_caucasian_sampleID.txt \
        --GMMATmodelFile=/mnt/plink/00_output/00_${1}_chrnc/out_${1}_chrnc.rda \
        --varianceRatioFile=/mnt/plink/00_output/00_${1}_chrnc/out_${1}_chrnc_cate.varianceRatio.txt \
        --SAIGEOutputFile=/mnt/plink/00_output/00_${1}_chrnc/assoc_${1}_chrnc_chr${2}.txt \
        --numLinesOutput=1 \
        --groupFile=/mnt/plink/01.efe_caucasian_LoF_MAF0.001_chr${2}.set.saige   \
        --sparseSigmaFile=/mnt/plink/00_output/00_${1}_chrnc/out_${1}_chrnc_cate.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE \
        --IsOutputAFinCaseCtrl=TRUE    \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsDropMissingDosages=TRUE \
        --LOCO=FALSE

# test_efe_caucasian_chr22_saige.set