## Alignment to Human Genome with TopHat2 ------------------------------------------------------------------------###

#!/bin/bash
#SBATCH --job-name=TopHat
#SBATCH --time=60:00:00
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --image=docker:araposo/tophat:1.0
#SBATCH --workdir=//mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/Reads/

export GENOME=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome
export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Tconv274 $GENOME AES-2/T274_Tconvs_conc1.fq.gz AES-2/T274_Tconvs_conc2.fq.gz

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Tconv276 $GENOME AES-4/T276_Tconvs_conc1.fq.gz AES-4/T276_Tconvs_conc2.fq.gz

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Tconv277 $GENOME AES-6/T277_Tconvs_conc1.fq.gz AES-6/T277_Tconvs_conc2.fq.gz

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Treg274 $GENOME AES-1/T274_Tregs_conc1.fq.gz AES-1/T274_Tregs_conc2.fq.gz

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Treg276 $GENOME AES-3/T276_Tregs_conc1.fq.gz AES-3/T276_Tregs_conc2.fq.gz

srun shifter tophat2 -g 1 -p 10  -G $ANNOTATION  -o TopHat_Treg277 $GENOME AES-5/T277_Tregs_conc1.fq.gz AES-5/T277_Tregs_conc2.fq.gz



echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID