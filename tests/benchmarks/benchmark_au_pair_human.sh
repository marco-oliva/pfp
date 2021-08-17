#!/bin/bash
#SBATCH --job-name=expAU-H
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --mem=1000gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_expAU-H.log
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
BASE_DIR_PD="${BASE_DIR_EXP}/human_tests"
PROFILER="/usr/bin/time --verbose"
AUPAIR="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/build/AuPair"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

d10_path="${BASE_DIR_PD}/30-07-2021_13-52-14"
d100_path="${BASE_DIR_PD}/30-07-2021_15-33-37"
d200_path="${BASE_DIR_PD}/31-07-2021_02-35-23"
d400_path="${BASE_DIR_PD}/01-08-2021_00-02-46"
d800_path="${BASE_DIR_PD}/03-08-2021_13-33-12"
d1600_path="${BASE_DIR_PD}/03-08-2021_21-35-37"
d2500_path="${BASE_DIR_PD}/04-08-2021_15-03-16"

##----------------------------------------------------------
# Run

${PROFILER} ${AUPAIR} -w 10 -i "${d10_path}/pfp" -o "${d10_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d100_path}/pfp" -o "${d100_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d200_path}/pfp" -o "${d200_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d400_path}/pfp" -o "${d400_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d800_path}/pfp" -o "${d800_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d1600_path}/pfp" -o "${d1600_path}/pfp_removed_ts"
${PROFILER} ${AUPAIR} -w 10 -i "${d2500_path}/pfp" -o "${d2500_path}/pfp_removed_ts"




