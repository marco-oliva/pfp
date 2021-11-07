#!/bin/bash
#SBATCH --job-name=PFPDCC
#SBATCH --account=boucher
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=1000gb
#SBATCH --time=240:00:00
#SBATCH --output=%j_PFP_FASTA_DCC.log
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

BENCHMARK_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/benchmarks_dcc22/benchmark_chr19.py"
DATA_HUMAN="/blue/boucher/marco.oliva/projects/experiments/pfp/DCC22/vcf_to_fa/data/samples"
SAMPLES_LIST="/blue/boucher/marco.oliva/projects/experiments/pfp/repo/pfp/tests/benchmarks/benchmarks_dcc22/input_list_2500.txt"

module load python/3.6
module load htslib
module load bcftools
module load git
module load gcc/9.3.0

${PROFILER} ${BENCHMARK_HUMAN} -m 1000 -d ${DATA_HUMAN} -s ${SAMPLES_LIST}
