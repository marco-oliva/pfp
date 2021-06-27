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
PROFILER="/usr/bin/time --verbose"
AUPAIR="/blue/boucher/marco.oliva/projects/experiments/pfp/AuPair/aupair"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

25_path="${BASE_DIR_EXP}/arabidopsis_tests/05-06-2021_09-09-31/pfp"
125_path="${BASE_DIR_EXP}/arabidopsis_tests/05-06-2021_10-30-33/pfp"
250_path="${BASE_DIR_EXP}/arabidopsis_tests/05-06-2021_13-01-29/pfp"
500_path="${BASE_DIR_EXP}/arabidopsis_tests/05-06-2021_16-59-18/pfp"
1000_path="${BASE_DIR_EXP}/arabidopsis_tests/06-06-2021_00-18-48/pfp"

##----------------------------------------------------------
# Run

${PROFILER} python ${AUPAIR} -w 10 -t 1000 - b 100 ${25_path}
#${PROFILER} python ${AUPAIR} -w 10 -t 1000 - b 100 ${125_path}
#${PROFILER} python ${AUPAIR} -w 10 -t 10000 - b 100 ${250_path}
#${PROFILER} python ${AUPAIR} -w 10 -t 10000 - b 100 ${500_path}
#${PROFILER} python ${AUPAIR} -w 10 -t 10000 - b 100 ${1000_path}



