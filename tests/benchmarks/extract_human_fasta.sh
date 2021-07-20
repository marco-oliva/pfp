#!/bin/bash
#SBATCH --job-name=ext-H
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_ext-H.log
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
PROFILER="/usr/bin/time --verbose"

EXTRACT_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/extract_human.py"
SAMPLES_LIST_BASE_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/samples_lists/input_list"
OUT_DIR="/blue/boucher/marco.oliva/projects/experiments/pfp/human_tests/data/samples"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

${PROFILER} ${EXTRACT_HUMAN} -t 32 -s "${SAMPLES_LIST_BASE_HUMAN}_1600.txt" -o ${OUT_DIR}
