l#!/bin/bash

#BSUB -P acc_kennylab
#BSUB -n 10
#BSUB -W 8:00
#BSUB -q premium
#BSUB -R span[hosts=1] 
#BSUB -R rusage[mem=2000]
#BSUB -J "gcta[1]"
#BSUB -o /sc/arion/projects/igh/kennylab/christa/heritability/logs/gcta.out.%J.%I 
#BSUB -e /sc/arion/projects/igh/kennylab/christa/heritability/logs/gcta.err.%J.%I 

. /sc/arion/projects/igh/kennylab/christa/miniforge/etc/profile.d/conda.sh
conda activate saige

chr=$LSB_JOBINDEX 

genotyping_data="../../data/biome/geno/gght_v2_topmed_allchr"
pop="european_all" 

../plink2 \
    --bfile $genotyping_data \
    --keep "fam/"$pop".txt" \
    --make-bed \
    --out "geno/"$pop 

gcta64 \
    --bfile "geno/"$pop  \
    --autosome \
    --maf 0.01 \
    --make-grm \
    --out "grm/"$pop  \
    --thread-num 10

# pheno="bmi"

# gcta64 \
#     --grm "grm/"$pop \
#     --pheno "pheno/"$pheno".txt" \
#     --reml \
#     --out "output/"$pop"_"$pheno \
#     --thread-num 10

# python calculate_average_cov.py