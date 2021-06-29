#!/bin/bash
#SBATCH --job-name=aupbA
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_aupbA.log
#SBATCH --constraint='hpg3&amd&milan&infiniband'
#
# Asking for hpg-milan 	64 	8 	8 	1 	512 	hpg3;amd;milan;infiniband 	AMD EPYC 75F3 32-Core Processor

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
BASE_DIR_PD="${BASE_DIR_EXP}/arabidopsis_tests"
PROFILER="/usr/bin/time --verbose"
AUPAIR="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/AuPair/aupair"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

d25_path="${BASE_DIR_PD}/05-06-2021_09-09-31"
d125_path="${BASE_DIR_PD}/05-06-2021_10-30-33"
d250_path="${BASE_DIR_PD}/05-06-2021_13-01-29"
d500_path="${BASE_DIR_PD}/05-06-2021_16-59-18"
d1000_path="${BASE_DIR_PD}/06-06-2021_00-18-48"

##----------------------------------------------------------
# Run

${PROFILER} python ${AUPAIR} -w 10 -t 100 -b 100 ${d25_path}/pfp
${PROFILER} python ${AUPAIR} -w 10 -t 100 -b 100 ${d125_path}/pfp
${PROFILER} python ${AUPAIR} -w 10 -t 100 -b 100 ${d250_path}/pfp
${PROFILER} python ${AUPAIR} -w 10 -t 1000 -b 100 ${d500_path}/pfp
${PROFILER} python ${AUPAIR} -w 10 -t 1000 -b 100 ${d1000_path}/pfp



