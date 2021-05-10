#!/bin/bash
#SBATCH --job-name=exfa
#SBATCH --account=boucher
#SBATCH --qos=boucher-b
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=%j_exfa.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
VCFTOFA="/blue/boucher/marco.oliva/projects/pfp/build/vcf_to_fa"
PROFILER="/usr/bin/time --verbose"

# Extract fasta, 1000 sequences
#${PROFILER} ${VCFTOFA} -v ${BASE_DIR_EXP}/data/vcf/chr19.vcf.gz -r ${BASE_DIR_EXP}/data/reference/19.fa.gz -m 1000 -o ${BASE_DIR_EXP}/data/fasta/19.1000.fa
#${PROFILER} ${VCFTOFA} -v ${BASE_DIR_EXP}/data/vcf/chr19.vcf.gz -r ${BASE_DIR_EXP}/data/reference/20.fa.gz -m 1000 -o ${BASE_DIR_EXP}/data/fasta/20.1000.fa
#${PROFILER} ${VCFTOFA} -v ${BASE_DIR_EXP}/data/vcf/chr19.vcf.gz -r ${BASE_DIR_EXP}/data/reference/21.fa.gz -m 1000 -o ${BASE_DIR_EXP}/data/fasta/21.1000.fa
${PROFILER} ${VCFTOFA} -v "${BASE_DIR_EXP}/data/vcf/chr19.vcf.gz" "${BASE_DIR_EXP}/data/vcf/chr20.vcf.gz" "${BASE_DIR_EXP}/data/vcf/chr21.vcf.gz" -r "${BASE_DIR_EXP}/data/reference/19.fa.gz" "${BASE_DIR_EXP}/data/reference/20.fa.gz" "${BASE_DIR_EXP}/data/reference/21.fa.gz" -m 1000 -o "${BASE_DIR_EXP}/data/fasta/19-21.1000.fa"



