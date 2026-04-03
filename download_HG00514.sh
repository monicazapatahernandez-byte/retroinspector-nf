#!/bin/bash
#SBATCH --job-name=dl_HG00514
#SBATCH --output=logs/dl_HG00514_%j.log
#SBATCH --error=logs/dl_HG00514_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

cd ~/tfm/simulacion/input_data

echo "=== Descargando HG00514 ==="
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00514.ONT.rebasecalled.fastq.tar.gz

echo "=== Extrayendo HG00514 ==="
tar -xzf HG00514.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > HG00514.fastq.gz
rm -rf projects/
rm -f HG00514.ONT.rebasecalled.fastq.tar.gz

echo "=== HG00514 listo ==="
ls -lh HG00514.fastq.gz
zcat HG00514.fastq.gz | head -4
du -sh ~/
