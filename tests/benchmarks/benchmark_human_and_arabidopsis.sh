#!/bin/bash
#SBATCH --job-name=ex-HeA
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_ex-HeA.log
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

BENCHMARK_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/benchmark_human.py"
SAMPLES_LIST_BASE_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/samples_lists/input_list"
DATA_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/human_tests/data/samples/samples"

BENCHMARK_ARABIDOPSIS="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/benchmark_arabidopsis.py"
SAMPLES_LIST_BASE_ARABIDOPSIS="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/samples_lists/arabidopsis_input_list"
DATA_ARABIDOPSIS="/blue/boucher/marco.oliva/projects/experiments/pfp/arabidopsis_tests/data/samples/samples"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_10.txt"
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_100.txt"
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_200.txt"
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_400.txt"
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_800.txt"
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_1600.txt" --skip-pscan
#${PROFILER} ${BENCHMARK_HUMAN} -t 32 -d ${DATA_HUMAN} -s "${SAMPLES_LIST_BASE_HUMAN}_all.txt" --skip-pscan

${PROFILER} ${BENCHMARK_ARABIDOPSIS} -t 32 -d ${DATA_ARABIDOPSIS} -s "${SAMPLES_LIST_BASE_ARABIDOPSIS}_25.txt"
${PROFILER} ${BENCHMARK_ARABIDOPSIS} -t 32 -d ${DATA_ARABIDOPSIS} -s "${SAMPLES_LIST_BASE_ARABIDOPSIS}_125.txt"
${PROFILER} ${BENCHMARK_ARABIDOPSIS} -t 32 -d ${DATA_ARABIDOPSIS} -s "${SAMPLES_LIST_BASE_ARABIDOPSIS}_250.txt"
${PROFILER} ${BENCHMARK_ARABIDOPSIS} -t 32 -d ${DATA_ARABIDOPSIS} -s "${SAMPLES_LIST_BASE_ARABIDOPSIS}_500.txt"
${PROFILER} ${BENCHMARK_ARABIDOPSIS} -t 32 -d ${DATA_ARABIDOPSIS} -s "${SAMPLES_LIST_BASE_ARABIDOPSIS}_1000.txt"

