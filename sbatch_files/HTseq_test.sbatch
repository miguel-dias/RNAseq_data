#!/bin/bash
#SBATCH --job-name=HT_Tc274
#SBATCH --time=7:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --workdir=//mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results/TopHat_Tconv274

export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf

##------- Sorts the Bam file ------------------------------------------------------------- ###

srun shifter --image=docker:argrosso/htstools samtools sort -@ 10 -n -o sortedbyname_TopHat_Tconv274.bam accepted_hits.bam

srun shifter --image=docker:argrosso/htstools samtools view -@ 10 -o sortedbyname_TopHat_Tconv274.sam  sortedbyname_TopHat_Tconv274.bam

srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 sortedbyname_TopHat_Tconv274.sam $ANNOTATION > TopHat_Tconv274.count

##srun shifter samtools view -o Tconv276.sam  accepted_hitsDT_sorted.bam

##srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 Tconv276.sam $ANNOTATION > Tconv276.count


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

