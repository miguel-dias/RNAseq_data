#!/bin/bash
#SBATCH --job-name=TopHat
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --nodes=2
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=8
#SBATCH --workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/

export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome
export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf
export Clean=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/PBMC_RNAseq_data/F17KTSEUHT0774_HUMjugE/Clean
export srun="srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog /mnt/beegfs/scratch/PRECISE/miguel.dias/$SLURM_JOB_ID.log"

export workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/
##-----------To run in selected folders----------##


ls -d $Clean/*/ | $parallel '$srun shifter --image=docker:araposo/tophat:1.0 tophat2 -g 1 -p $SLURM_CPUS_PER_TASK -G $ANNOTATION -o $workdir/$(basename {}) $GENOME {}*_conc1.fq.gz {}*_conc2.fq.gz'


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
