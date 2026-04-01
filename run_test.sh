#!/bin/bash
#SBATCH --job-name=retroinspector-nf
#SBATCH --output=logs/retroinspector_%j.log
#SBATCH --error=logs/retroinspector_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/alumno27/miniconda3/envs/nf-core

cd ~/tfm/retroinspector-nf

nextflow run main.nf \
    --input test_data/samplesheet.csv \
    --reference /home/software/nanopore/hg38.fa \
    --outdir results_test \
    --threads 4 \
    --mode full \
    -profile dayhoff \
    -resume
