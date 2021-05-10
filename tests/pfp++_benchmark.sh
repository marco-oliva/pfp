#!/bin/bash
#SBATCH --job-name=timep
#SBATCH --account=boucher
#SBATCH --qos=boucher-b
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=%j_timep.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
PFP="/blue/boucher/marco.oliva/projects/experiments/pfp/tmp/bin/pfp++"
NEWSCAN="/blue/boucher/marco.oliva/projects/experiments/pfp/tmp/bin/newscan.x"
PROFILER="/usr/bin/time --verbose"

# Extract fasta, 1000 sequences
#${PROFILER} ${PFP} -w 10 -p 100 -t 32 -m 1000 --use-acceleration -v ${BASE_DIR_EXP}/data/vcf/chr19.vcf.gz -r ${BASE_DIR_EXP}/data/reference/19.fa.gz -o ${BASE_DIR_EXP}/tmp/19.out
#${PROFILER} ${PFP} -w 10 -p 100 -t 32 -m 1000 --use-acceleration -v ${BASE_DIR_EXP}/data/vcf/chr20.vcf.gz -r ${BASE_DIR_EXP}/data/reference/20.fa.gz -o ${BASE_DIR_EXP}/tmp/20.out
#${PROFILER} ${PFP} -w 10 -p 100 -t 32 -m 1000 --use-acceleration -v ${BASE_DIR_EXP}/data/vcf/chr21.vcf.gz -r ${BASE_DIR_EXP}/data/reference/21.fa.gz -o ${BASE_DIR_EXP}/tmp/21.out


${PROFILER} ${PFP} --configure ${BASE_DIR_EXP}/tmp/config.ini -w 10 -p 100 -t 32 --use-acceleration -m 1000 -o ${BASE_DIR_EXP}/tmp/19-21.out

${PROFILER} ${NEWSCAN} ${BASE_DIR_EXP}/data/fasta/19-21.1000.fa -f -w 10 -p 100 -t 32


