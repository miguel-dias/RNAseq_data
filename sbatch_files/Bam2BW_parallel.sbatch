#!/bin/bash
#SBATCH --job-name=Bam2BW
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=6
#SBATCH --image=argrosso/htstools:latest
#SBATCH --workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/Bam_2_BigWig/

export TopHat_results=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results
export srun="srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog /mnt/beegfs/scratch/PRECISE/miguel.dias/Bam_2_BigWig/$SLURM_JOB_ID.log"

export workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/Bam_2_BigWig

## --------- Sort the bam file, creates index and converts to bigwig for USCS Genome Browser----------------------------------------- ###

ls -d $TopHat_results/*/ | $parallel '$srun shifter samtools sort -@ $SLURM_CPUS_PER_TASK -o $workdir/sorted_$(basename {}).bam {}accepted_hits.bam; $srun shifter samtools index $workdir/sorted_$(basename {}).bam $workdir/sorted_$(basename {}).bam.bai; $srun shifter bamCoverage -b $workdir/sorted_$(basename {}).bam -o $workdir/sorted_$(basename {}).bw -of bigwig --binSize 1 --numberOfProcessors max/2'

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

