#!/bin/bash -l

#SBATCH --job-name=subset_vcf_pop
#SBATCH -o subset_vcf_pop.%A.%a.out
#SBATCH -e subset_vcf_pop.%A.%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --array=1-7
#SBATCH --mail-user=melodysyue@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --chdir=/scratch-lustre/peuclide/Sockeye/raw_data_1/yue

enable_lmod
module load vcftools

POP=$(find ./popmap -name '*.list' | sort |awk " NR==$SLURM_ARRAY_TASK_ID")
NAME=$(basename ${POP} .list)

echo ${POP}

vcftools --gzvcf LG1-29_output.recalibrated.filtered.1.no-indels.no-fail.biallele.miss.0.9.recode.vcf.gz \
	--recode --recode-INFO-all --keep $POP --out ./vcf_pop/$NAME
