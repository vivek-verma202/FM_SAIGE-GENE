# shellcheck shell=bash
mkdir /home/vverma3/repo/FM_SAIGE-GENE/data

# only chop SNPs, not samples
for chr in {1..22}; do
/home/vverma3/plink2 \
  --threads 99   \
  --bfile /mnt/nfs/backup/data/uk_biobank/"${chr}" \
  --maf  0.05    \
  --geno 0.05    \
  --hwe 1e-12    \
  --make-bed     \
  --out /home/vverma3/repo/FM_SAIGE-GENE/data/mini"${chr}"
done

# Then use plink again to merge all the minis:
# create list of files to merge
for chr in {1..22}; do
  echo "/home/vverma3/repo/FM_SAIGE-GENE/data/mini""${chr}"
done > merge.list

# merge the files
/home/mparis15/PLINK/plink \
  --threads 99 \
  --merge-list merge.list \
  --make-bed \
  --out /home/vverma3/repo/FM_SAIGE-GENE/data/merged

# clean
rm -r /home/vverma3/repo/FM_SAIGE-GENE/data/*
rm /home/vverma3/repo/FM_SAIGE-GENE/merge.list \
/home/vverma3/repo/FM_SAIGE-GENE/merged.log \
/home/vverma3/repo/FM_SAIGE-GENE/merged.nosex \
/home/vverma3/repo/FM_SAIGE-GENE/plink2.log
mv /home/vverma3/repo/FM_SAIGE-GENE/merged.* /home/vverma3/repo/FM_SAIGE-GENE/data

# create GRM
/opt/R/3.6.0/bin/Rscript ./SAIGE/extdata/createSparseGRM.R \
--plinkFile=/home/vverma3/repo/FM_SAIGE-GENE/data/merged   \
--nThreads=108                                             \
--outputPrefix=/home/vverma3/repo/FM_SAIGE-GENE/data/GRM   \
--minMAFforGRM=0.05


