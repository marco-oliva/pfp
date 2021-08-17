#!/bin/bash
#SBATCH --job-name=expAU-A
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=1000gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_expAU-A.log
#SBATCH --constraint='hpg3&amd&rome&infiniband'
#
# Asking for hpg-default 	128 	8 	16 	1 	1028 	hpg3;amd;rome;infiniband 	AMD EPYC 7702 64-Core Processor

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date;

if command -v lshw &> /dev/null
then
  lshw
else
  echo "lshw could not be found"
fi

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
BASE_DIR_PD="${BASE_DIR_EXP}/arabidopsis_tests"
PROFILER="/usr/bin/time --verbose"
AUPAIR="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/build/AuPair"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

d25_path="${BASE_DIR_PD}/28-07-2021_09-32-03"
d125_path="${BASE_DIR_PD}/28-07-2021_09-43-23"
d250_path="${BASE_DIR_PD}/28-07-2021_10-20-16"
d500_path="${BASE_DIR_PD}/28-07-2021_11-29-09"
d1000_path="${BASE_DIR_PD}/28-07-2021_13-47-20"

##----------------------------------------------------------
# Run

${PROFILER} ${AUPAIR} -w 10 -i "${d25_path}/pfp" -o "${d25_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d125_path}/pfp" -o "${d125_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d250_path}/pfp" -o "${d250_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d500_path}/pfp" -o "${d500_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d1000_path}/pfp" -o "${d1000_path}/pfp_removed_ts"





