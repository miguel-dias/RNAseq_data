#!/bin/bash
#SBATCH --job-name=FW_test
#SBATCH --time=2:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=argrosso/htstools:latest
#SBATCH --workdir=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results/TopHat_Tconv274

## FW

srun shifter bamCoverage -b accepted_hitsDT_sorted_position.bam -o testFW.bw --filterRNAstrand forward
srun shifter bamCoverage -b accepted_hitsDT_sorted_position.bam -o testRV.bw --filterRNAstrand reverse
