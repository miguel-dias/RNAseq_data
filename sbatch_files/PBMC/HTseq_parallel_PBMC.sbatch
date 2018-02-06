#!/bin/bash
#SBATCH --job-name=HTpara
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --nodes=3
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=12
#SBATCH --workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/PBMC_HTseq_results/

export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf
export PBMC_TopHat_results=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/PBMC_TopHat_results

export srun="srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog /mnt/beegfs/scratch/PRECISE/miguel.dias/results_HTseq/$SLURM_JOB_ID.log"

##para ver dentro das plicas, faz-se o export da wk##
export workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/PBMC_HTseq_results

ls -d $PBMC_TopHat_results/*/ | $parallel '$srun shifter --image=docker:argrosso/htstools samtools sort -@ $SLURM_CPUS_PER_TASK -n -o $workdir/sortedbyname_$(basename {}).bam {}accepted_hits.bam; $srun shifter --image=docker:argrosso/htstools samtools view -@ $SLURM_CPUS_PER_TASK -o $workdir/sortedbyname_$(basename {}).sam  $workdir/sortedbyname_$(basename {}).bam; $srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 $workdir/sortedbyname_$(basename {}).sam $ANNOTATION > $workdir/$(basename {}).count' 

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID