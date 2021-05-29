#!/bin/bash
#SBATCH --job-name=exhbp
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60gb
#SBATCH --time=72:00:00
#SBATCH --output=%j_exhbp.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
PROFILER="/usr/bin/time --verbose"
BENCHMARK="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/benchmark_human.py"
SAMPLES_LIST_BASE="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/input_list"

module load python
module load htslib
module load bcftools
module load git

${PROFILER} ${BENCHMARK} -t 32 -s "${SAMPLES_LIST_BASE}_10.txt"
${PROFILER} ${BENCHMARK} -t 32 -s "${SAMPLES_LIST_BASE}_100.txt"
${PROFILER} ${BENCHMARK} -t 32 -s "${SAMPLES_LIST_BASE}_200.txt"

