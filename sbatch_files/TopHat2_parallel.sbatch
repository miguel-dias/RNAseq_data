#!/bin/bash
#SBATCH --job-name=TopHat
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=6
#SBATCH --image=docker:araposo/tophat:1.0
#SBATCH --workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/Tophat/

export Thymus_RNAseq_data=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/Thymus_RNAseq_data/
export srun="srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
parallel="parallel --dry-run --delay 0.2 -j $SLURM_NTASKS --joblog /mnt/beegfs/scratch/PRECISE/miguel.dias/Tophat/$SLURM_JOB_ID.log"

export workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/Tophat

##-----------To run in selected folders----------##

export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome
export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf

$parallel 'srun shifter tophat2 -g 1 -p 6  -G $ANNOTATION  -o $workdir/$basename{}) $GENOME *_conc1.fq.gz *_conc2.fq.gz' ::: AES-1 AES-2 AES-3 AES-4 AES-5 AES-6
#$parallel 'srun shifter tophat2 -g 1 -p 6  -G $ANNOTATION  -o $workdir/$basename{}) $GENOME {#}_conc1.fq.gz {#}_conc2.fq.gz' ::: AES-1 AES-2 AES-3 AES-4 AES-5 AES-6

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
