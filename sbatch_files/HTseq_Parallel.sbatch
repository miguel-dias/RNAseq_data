#!/bin/bash
#SBATCH --job-name=para_test
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=10
#SBATCH --workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/results_HTseq/

export ANNOTATION=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Annotation/gencode.v26.annotation.gtf
export TopHat_results=/mnt/nfs/lobo/PRECISE-NFS/miguel.dias/Human_Genome/TopHat_results

## export wk=$SLURM_SUBMIT_DIR ## working directory


## export log_dir=$wk/TopHat_Tconv277
## export log_dir=$wk/TopHat_Tconv274
## export log_dir=$wk/TopHat_Treg274
## export log_dir=$wk/TopHat_Treg277

export srun="srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog /mnt/beegfs/scratch/PRECISE/miguel.dias/results_HTseq/$SLURM_JOB_ID.log"

##para ver dentro das plicas, faz-se o export da wk##
export workdir=/mnt/beegfs/scratch/PRECISE/miguel.dias/results_HTseq

ls -d $TopHat_results/*/ | $parallel '$srun shifter --image=docker:argrosso/htstools samtools sort -@ $SLURM_CPUS_PER_TASK -n -o $workdir/sortedbyname_$(basename {}).bam {}accepted_hits.bam; $srun shifter --image=docker:argrosso/htstools samtools view -@ $SLURM_CPUS_PER_TASK -o $workdir/sortedbyname_$(basename {}).sam  $workdir/sortedbyname_$(basename {}).bam; $srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 $workdir/sortedbyname_$(basename {}).sam $ANNOTATION > $workdir/$(basename {}).count' 

##-----------To run in selected folders----------##

##$parallel 'srun shifter --image=docker:argrosso/htstools samtools sort -n -o sortedbyname_{#}.bam $TopHat_results/${folder}/accepted_hits.bam; srun shifter --image=docker:argrosso/htstools samtools view -o sortedbyname.sam  sortedbyname.bam; srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 sortedbyname.sam $ANNOTATION > {/.}.count' ::: TopHat_Tconv277 TopHat_Tconv274 TopHat_Treg274 TopHat_Treg277


#for folder in TopHat_Tconv277 TopHat_Tconv274 TopHat_Treg274 TopHat_Treg277; do 
#
#  export folder
#  echo $folder
#  ls $TopHat_results/${folder}
#  time ls $TopHat_results/${folder} | $parallel 'srun shifter --image=docker:argrosso/htstools samtools sort -n -o sortedbyname_{#}.bam $TopHat_results/${folder}/accepted_hits.bam; srun shifter --image=docker:argrosso/htstools samtools view -o sortedbyname.sam  sortedbyname.bam; srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 sortedbyname.sam $ANNOTATION > {/.}.count'
#  
#  ## time ls $TopHat_results/${folder} | $parallel 'srun shifter --image=docker:argrosso/htstools samtools view -o sortedbyname.sam  sortedbyname.bam'
#  ## time ls $TopHat_results/${folder} | $parallel 'srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 sortedbyname.sam $ANNOTATION > {/.}.count'
#
#done

## ls ${folder}/*.bam | $parallel '$srun shifter samtools sort -o {/.}.sorted.bam {}; $srun shifter samtools index {}' 
## ls $log_dir/TopHat_Treg277 | $parallel 'srun shifter --image=docker:argrosso/htstools samtools sort -n -o sortedbyname.bam accepted_hits.bam'

## ls $log_dir/TopHat_Treg277 | $parallel 'srun shifter --image=docker:argrosso/htstools samtools view -o Treg277_name.sam  sortedbyname.bam'

## ls $log_dir/TopHat_Treg277 | $parallel 'srun shifter --image=araposo/tophat:1.0 htseq-count -s no -a 10 Treg277_name.sam $ANNOTATION > Treg277_name.count'


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID