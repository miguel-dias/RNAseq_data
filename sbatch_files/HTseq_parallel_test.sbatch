#!/bin/bash
#SBATCH --job-name=HT_test
#SBATCH --time=20:00:00
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:argrosso/htstools
#SBATCH --workdir=//mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results/

export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf

export wk=$SLURM_SUBMIT_DIR ## working directory

export log_dir=$wk/TopHat_Tconv276
export log_dir=$wk/TopHat_Tconv277
export log_dir=$wk/TopHat_Tconv274
export log_dir=$wk/TopHat_Treg276
export log_dir=$wk/TopHat_Treg274
export log_dir=$wk/TopHat_Treg277

srun="srun -N1 -n1"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/HTseq_count/$SLURM_JOB_ID.log"

for folder in TopHat_Tconv276 TopHat_Tconv277 TopHat_Tconv274
ls ${folder}/*.bam | $parallel '$srun shifter  samtools sort -o {/.}.sorted.bam {}; $srun shifter samtools index {}' 

##time ls $log_dir/TopHat_Tconv276 | $parallel 'srun shifter samtools view -o Tconv274.sam accepted_hitsDT_sorted.bam'

##time ls $log_dir/TopHat_Tconv276 | $parallel 'srun shifter htseq-count -s no -a 10 Tconv274.sam $ANNOTATION > Tconv274.count'


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
