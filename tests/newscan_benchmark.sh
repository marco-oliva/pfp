#!/bin/bash
#SBATCH --job-name=timen
#SBATCH --account=boucher
#SBATCH --qos=boucher-b
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=%j_timen.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date

##----------------------------------------------------------
# Setup
BASE_DIR_EXP="/blue/boucher/marco.oliva/projects/experiments/pfp"
NEWSCAN="/blue/boucher/marco.oliva/projects/experiments/pfp/tmp/bin/newscan.x"
PROFILER="/usr/bin/time --verbose"

# Extract fasta, 1000 sequences
${PROFILER} ${NEWSCAN} ${BASE_DIR_EXP}/data/fasta/19.1000.fa -f -w 10 -p 100 -t 32
${PROFILER} ${NEWSCAN} ${BASE_DIR_EXP}/data/fasta/20.1000.fa -f -w 10 -p 100 -t 32
${PROFILER} ${NEWSCAN} ${BASE_DIR_EXP}/data/fasta/21.1000.fa -f -w 10 -p 100 -t 32

