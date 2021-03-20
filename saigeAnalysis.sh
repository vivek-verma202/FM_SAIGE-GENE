# apply saige-gene to do rare variants analysis

cd /scratch/xao1/project/Rare_variant_pain/pheno

# remove samples with PHENO == 8 cause they are samples having general pain
awk '{if($NF != 8) print $0}' covars_and_pheno.txt > covars_and_pheno_fSITE_chrnc.txt

cd /scratch/xao1/project/Rare_variant_pain/saige_50k

awk '{print $1"\t"$1}' exclude_kinship.txt > 01.exclude_kinship_sampleID.txt

cat 01.exclude_kinship_sampleID.txt ukbb_withdrewSampleID_20210201.txt | sort | uniq > 01.kinship_retired_sampleID.txt

## prepare data

# prepare caucasian plink file

######### keep LoF variants

plink --bfile efe --keep caucasion_list.txt --remove ukbb_withdrewSampleID_20210201.txt --out 01.efe_caucasian --make-bed

plink --bfile 01.efe_caucasian --remove 01.kinship_retired_sampleID.txt --out 01.efe_caucasian_nokin --make-bed

plink --bfile 01.efe_caucasian_nokin --extract 01.LoF.MAF_0.001.variantID.txt --out 01.efe_caucasian_nokin --freq
# PLINK v1.90p 64-bit (28 Oct 2019)              www.cog-genomics.org/plink/1.9/
# (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to 01.efe_caucasian_nokin.log.
# Options in effect:
#   --bfile 01.efe_caucasian_nokin
#   --extract 01.LoF.MAF_0.001.variantID.txt
#   --freq
#   --out 01.efe_caucasian_nokin

# 257858 MB RAM detected; reserving 128929 MB for main workspace.
# 8959608 variants loaded from .bim file.
# 39009 people (18041 males, 20968 females) loaded from .fam.
# --extract: 73448 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 39009 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 5423 het. haploid genotypes present (see 01.efe_caucasian_nokin.hh );
# many commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.993719.
# --freq: Allele frequencies (founders only) written to
# 01.efe_caucasian_nokin.frq .

# prepare SNP list having MAF < 0.001 in unrelated ukb white population and form them in SAIGE format
cat 01.efe_caucasian_nokin.frq | sed 's/:/\t/g' | awk 'NR>1{if($8 < 0.001) print $1":"$3"_"$7"/"$6}' > 01.LoFs_MAF0.001_SAIGE_variantID.txt

# no duplicate variantID found
cat 01.LoFs_MAF0.001_SAIGE_variantID.txt | sort | uniq -c | awk '{if($1 > 1) print}' | wc -l 
# 0

# new file to update snp name: column1: old name, column2: new name(chr:pos_ref/alt)
cat 01.efe_caucasian.bim | awk -F "\t" '{print $2"\t"$1":"$4"_"$6"/"$5}' > 01.efe_caucasian.newbim
# found duplicated new names
cat 01.efe_caucasian.newbim | awk '{dup[$2]++; if(dup[$2] > 1) print $0"_DUP"dup[$2]; else print $0 }' > 01.efe_caucasian_uniq.newbim
# update snp name and make gene set 
plink --bfile 01.efe_caucasian --update-name 01.efe_caucasian_uniq.newbim --extract 01.LoFs_MAF0.001_SAIGE_variantID.txt --make-set gencode.v26.gene.list --write-set --out 01.efe_caucasian_LoF_MAF0.001
# PLINK v1.90p 64-bit (28 Oct 2019)              www.cog-genomics.org/plink/1.9/
# (C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to 01.efe_caucasian_LoF_MAF0.001.log.
# Options in effect:
#   --bfile 01.efe_caucasian
#   --extract 01.LoFs_MAF0.001_SAIGE_variantID.txt
#   --make-set gencode.v26.gene.list
#   --out 01.efe_caucasian_LoF_MAF0.001
#   --update-name 01.efe_caucasian_uniq.newbim
#   --write-set

# 257858 MB RAM detected; reserving 128929 MB for main workspace.
# Warning: Unusually long new variant ID(s) in --update-name file.  Double-check
# your file and command-line parameters, and consider changing your naming
# scheme if you encounter memory problems.
# 8959608 variants loaded from .bim file.
# 41244 people (18971 males, 22273 females) loaded from .fam.
# --update-name: 8959608 values updated.
# --extract: 72613 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 41244 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Warning: 2705 het. haploid genotypes present (see
# 01.efe_caucasian_LoF_MAF0.001.hh ); many commands treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
# Total genotyping rate is 0.994579.
# --make-set: 58219 sets defined.
# 72613 variants and 41244 people pass filters and QC.
# Note: No phenotypes present.
# --write-set: 01.efe_caucasian_LoF_MAF0.001.set written.

# make new set as saige step2 input
cat 01.efe_caucasian_LoF_MAF0.001.set | sed '/DUP/d' | sed 's/END//' | sed '/^$/d' | tr '\n' '\t' | sed 's/ENSG/\nENSG/g' | sed '/^$/d; s/\t$//g' | awk -F "\t" '{if(NF>2) print}' > 01.efe_caucasian_LoF_MAF0.001.set.saige

# make set for each chromosome
for chr in `seq 1 22`; do
  plink --bfile 01.efe_caucasian --chr $chr --update-name 01.efe_caucasian_uniq.newbim --extract 01.LoFs_MAF0.001_SAIGE_variantID.txt --make-set gencode.v26.gene.list --write-set --out 01.efe_caucasian_LoF_MAF0.001_chr${chr}
  cat 01.efe_caucasian_LoF_MAF0.001_chr${chr}.set | sed '/DUP/d' | sed 's/END//' | sed '/^$/d' | tr '\n' '\t' | sed 's/ENSG/\nENSG/g' | sed '/^$/d; s/\t$//g' | awk -F "\t" '{if(NF>2) print}' > 01.efe_caucasian_LoF_MAF0.001_chr${chr}.set.saige
done

rm 01.efe_caucasian_LoF_MAF0.001_chr[1-9].set 01.efe_caucasian_LoF_MAF0.001_chr[1-9].log 
rm 01.efe_caucasian_LoF_MAF0.001_chr[1-2][0-9].set 01.efe_caucasian_LoF_MAF0.001_chr[1-2][0-9].log

# make sampleID file as saige step2 input
cat 01.efe_caucasian.fam | cut -d ' ' -f1 > 01.efe_caucasian_sampleID.txt

# convert plink format to vcf format as saige step2 input
seq 24 -1 1 | while read chrom; do
	plink --bfile 01.efe_caucasian --chr $chrom --update-name 01.efe_caucasian_uniq.newbim  --recode vcf-iid bgz --out 01.efe_caucasian_${chrom}
	tabix -f -p vcf 01.efe_caucasian_${chrom}.vcf.gz
done 

##############################################################
# do step0 : create sparse GRM
nice -n 10 Rscript /usr/local/bin/createSparseGRM.R \
	--plinkFile=01.efe_caucasian \
	--nThreads=28  \
	--outputPrefix=./00_output/efe_caucasian_sparseGRM	\
	--numRandomMarkerforSparseKin=2000	\
	--relatednessCutoff=0.125 

##############################################################
# do step1 for chronic pain 

pain=chrnc
for site in GeneralPain NeckShoulder HipPain BackPain StomachAbdominal KneePain Headache FacialPain; do
    cd /scratch/xao1/project/Rare_variant_pain/saige_50k
    if [ ! -d ./00_output/00_${site}_${pain} ]; then
    	mkdir ./00_output/00_${site}_${pain}
    fi


#---------- build the covars+pheno file
echo "FID,IID,SEX,AGE,AGE_SQ,UKBL,PCA1,PCA2,PCA3,PCA4,PCA5,PCA6,PCA7,PCA8,PCA9,PCA10,PCA11,PCA12,PCA13,PCA14,PCA15,PCA16,PCA17,PCA18,PCA19,PCA20,PCA21,PCA22,PCA23,PCA24,PCA25,PCA26,PCA27,PCA28,PCA29,PCA30,PCA31,PCA32,PCA33,PCA34,PCA35,PCA36,PCA37,PCA38,PCA39,PCA40,C_10003,C_11001,C_11002,C_11003,C_11004,C_11005,C_11006,C_11007,C_11008,C_11009,C_11011,C_11012,C_11013,C_11014,C_11016,C_11017,C_11018,C_11020,C_11021,C_11022,C_11023,PHENO" \
  | tr ',' $'\t' \
  > /scratch/xao1/project/Rare_variant_pain/pheno/covars_and_pheno_${site}_${pain}.txt
join -t $'\t' -1 1 -2 1 \
<( \
  tail -n +2 /scratch/xao1/project/Rare_variant_pain/pheno/covars_full.txt \
  | LANG=en_EN sort -t $'\t' -k 1,1 \
 ) \
<( \
  tail -n +2 /scratch/xao1/project/Rare_variant_pain/pheno/pheno_${site}_${pain}.txt \
  | cut -d $'\t' -f 1,3 \
  | LANG=en_EN sort -t $'\t' -k 1,1 \
 ) \
 >> /scratch/xao1/project/Rare_variant_pain/pheno/covars_and_pheno_${site}_${pain}.txt

if [ "${site}" == "SITE" ]; then
time \
  /home/xao1/project/Rare_variant_pain/script/run_saige_step1.sh \
  quantitative $site
else
time \
  /home/xao1/project/Rare_variant_pain/script/run_saige_step1.sh \
  binary $site
fi
done

# do step 1 for multisite pain   
docker run --user $(id -u):$(getent group ukb | cut -d ':' -f 3) \
     	-v /scratch/xao1/project/Rare_variant_pain/saige_50k:/mnt/plink \
     	-v /scratch/xao1/project/Rare_variant_pain:/mnt/ \
     	-w /scratch/xao1/project/Rare_variant_pain/saige_50k \
		6eef0eeb60e2 nice -n 10 Rscript /usr/local/bin/step1_fitNULLGLMM.R     \
		--plinkFile=/mnt/plink/01.efe_caucasian \
        --phenoFile=/mnt/pheno/covars_and_pheno_fSITE_chrnc.txt \
        --phenoCol=PHENO \
        --covarColList=SEX,AGE,AGE_SQ,UKBL,PCA1,PCA2,PCA3,PCA4,PCA5,PCA6,PCA7,PCA8,PCA9,PCA10,PCA11,PCA12,PCA13,PCA14,PCA15,PCA16,PCA17,PCA18,PCA19,PCA20,PCA21,PCA22,PCA23,PCA24,PCA25,PCA26,PCA27,PCA28,PCA29,PCA30,PCA31,PCA32,PCA33,PCA34,PCA35,PCA36,PCA37,PCA38,PCA39,PCA40,C_10003,C_11001,C_11002,C_11003,C_11004,C_11005,C_11006,C_11007,C_11008,C_11009,C_11011,C_11012,C_11013,C_11014,C_11016,C_11017,C_11018,C_11020,C_11021,C_11022,C_11023 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
        --invNormalize=TRUE     \
        --outputPrefix=/mnt/plink/00_output/00_multisite_pain/out_multisite \
		--outputPrefix_varRatio=/mnt/plink/00_output/00_multisite_pain/out_multisite_cate \
		--sparseGRMFile=/mnt/plink/00_output/efe_caucasian_sparseGRM_numRandomMarkerforSparseKin_2000_relatednessCutoff_0.125_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx    \
        --sparseGRMSampleIDFile=/mnt/plink/00_output/efe_caucasian_sparseGRM_numRandomMarkerforSparseKin_2000_relatednessCutoff_0.125_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
        --nThreads=28 \
        --LOCO=FALSE \
		--skipModelFitting=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE	

##############################################################
# do step2 

# A test
plink --bfile 01.efe_caucasian --chr 22 --update-name 01.efe_caucasian_uniq.newbim  --make-set gencode.v26.gene.list --write-set --out test_efe_caucasian_chr22

cat test_efe_caucasian_chr22.set | sed '/DUP/d' | sed 's/END//' | sed '/^$/d' | tr '\n' '\t' | sed 's/ENSG/\nENSG/g' | sed '/^$/d' | awk -F "\t" '{if(NF>2) print}' > test_efe_caucasian_chr22_saige.set

pain=chrnc
site=multisite
chr=22
cd /scratch/xao1/project/Rare_variant_pain/saige_50k
time /home/xao1/project/Rare_variant_pain/script/run_saige_step2.sh $site $chr

pain=chrnc
for site in multisite GeneralPain NeckShoulder HipPain BackPain StomachAbdominal KneePain Headache FacialPain; do
    cd /scratch/xao1/project/Rare_variant_pain/saige_50k
    
    for chr in `seq 22 -1 1`; do
    	time /home/xao1/project/Rare_variant_pain/script/run_saige_step2.sh $site $chr
	    done
done


