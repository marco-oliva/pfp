#!/bin/bash
#SBATCH --job-name=pdisp
#SBATCH --account=boucher
#SBATCH --qos=boucher
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=marco.oliva@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=96:00:00
#SBATCH --output=logs/%j_disp.log
#SBATCH --error=logs/%j_disp.log

##----------------------------------------------------------
# Print Some statistics
pwd; hostname; date


##----------------------------------------------------------
# Modules
module load snakemake


##----------------------------------------------------------
# Run

snakemake --cluster "sbatch -A {cluster.account} -q {cluster.qos} -c {cluster.cpus-per-task} -N {cluster.Nodes}  -t {cluster.runtime} --mem {cluster.mem} -J {cluster.jobname} --mail-type={cluster.mail_type} --mail-user={cluster.mail} --output {cluster.out} --error {cluster.err}" --cluster-config cluster.json --jobs 100 --latency-wait 20 --rerun-incomplete
