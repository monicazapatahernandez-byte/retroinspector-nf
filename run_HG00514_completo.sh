#!/bin/bash
#SBATCH --job-name=retro_HG00514
#SBATCH --output=logs/retro_HG00514_%j.log
#SBATCH --error=logs/retro_HG00514_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/alumno27/miniconda3/envs/nf-core

DATA_DIR=~/tfm/retroinspector-nf/input_data
PIPE_DIR=~/tfm/retroinspector-nf

echo "=== INICIO: $(date) ==="
du -sh ~

# Limpiar posibles restos
rm -f $DATA_DIR/HG00514.ONT.rebasecalled.fastq.tar.gz
rm -rf $DATA_DIR/projects/
rm -f $DATA_DIR/HG00514.fastq.gz

echo "=== Descargando HG00514 ==="
cd $DATA_DIR
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00514.ONT.rebasecalled.fastq.tar.gz

echo "=== Verificando integridad ==="
if ! gzip -t HG00514.ONT.rebasecalled.fastq.tar.gz 2>/dev/null; then
    echo "ERROR: tar.gz corrupto"
    exit 1
fi

echo "=== Extrayendo ==="
tar -xzf HG00514.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > HG00514.fastq.gz
rm -rf projects/
rm -f HG00514.ONT.rebasecalled.fastq.tar.gz

echo "=== FASTQ listo: $(ls -lh HG00514.fastq.gz) ==="
zcat HG00514.fastq.gz | head -4
du -sh ~

echo "=== Lanzando Nextflow ==="
cd $PIPE_DIR
/home/alumno27/miniconda3/envs/nf-core/bin/nextflow run main.nf \
    --input test_data/samplesheet_HG00514.csv \
    --reference /home/software/nanopore/hg38.fa \
    --outdir results_HG00514 \
    -profile dayhoff \
    -resume

echo "=== FIN: $(date) ==="
du -sh ~
