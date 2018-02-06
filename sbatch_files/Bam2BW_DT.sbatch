#!/bin/bash
#SBATCH --job-name=FwdRev_test
#SBATCH --time=30:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=argrosso/htstools:latest
#SBATCH --workdir=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results/TopHat_Tconv277


## ------ To get the file for transcripts that originated from the forward strand -------##
 
##------- Sorts the Bam file ------------------------------------------------------------- ###

srun shifter samtools sort -o accepted_hitsDT_sorted.bam accepted_hits.bam 

## ------- Include reads that are 2nd in a pair (128) and exclude reads that are mapped to the reverse strand (16)------ ##

srun shifter samtools view -b -f 128 -F 16 accepted_hitsDT_sorted.bam > accepted_hitsDT_sorted.fwd1.bam

## --------- Exclude reads that are mapped to the reverse strand (16) and first in a pair (64): 64+16=80 --------- ##

srun shifter samtools view -b -f 80 accepted_hitsDT_sorted.bam > accepted_hitsDT_sorted.fwd2.bam

##---------- Combine the temporary files---------##
srun shifter samtools merge -f fwd.bam accepted_hitsDT_sorted.fwd1.bam accepted_hitsDT_sorted.fwd2.bam 

## ------ index the filtered BAM file ------##

srun shifter samtools index fwd.bam 

##-------- run bam Coverage --------##

srun shifter bamCoverage -b fwd.bam -o accepted_hitsDT_sorted.fwd.bigWig -of bigwig --binSize 1 --numberOfProcessors max/2

## ----remove temporary files-----##
srun shifter rm accepted_hitsDT_sorted.fwd*.bam

## ------ To get the file for transcripts that originated from the reverse strand -------##

##---include reads that map to the reverse strand (128) and are second in a pair (16): 128+16=144)----##

srun shifter samtools view -b -f 144 accepted_hitsDT_sorted.bam > accepted_hitsDT_sorted.rev1.bam

##----include reads that are first in a pair (64) but exclude those ones that map to the reverse strand (16)----##
srun shifter samtools view -b -f 64 -F 16 accepted_hitsDT_sorted.bam > accepted_hitsDT_sorted.rev2.bam

##---Merge the temporary files ---##
srun shifter samtools merge -f rev.bam accepted_hitsDT_sorted.rev1.bam accepted_hitsDT_sorted.rev2.bam

##----index the merged, filtered BAM file----##
srun shifter samtools index rev.bam

##----run bamCoverage-----##
srun shifter bamCoverage -b rev.bam -o accepted_hitsDT_sorted.rev.bw -of bigwig --binSize 1 --numberOfProcessors max/2

##----remove temporary files----##
srun shifter rm accepted_hitsDT_sorted.rev*.bam


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
