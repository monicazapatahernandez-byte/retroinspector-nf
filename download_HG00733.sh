#!/bin/bash
#SBATCH --job-name=dl_HG00733
#SBATCH --output=logs/dl_HG00733_%j.log
#SBATCH --error=logs/dl_HG00733_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

cd ~/tfm/simulacion/input_data

echo "=== Descargando HG00733 ==="
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00733.ONT.rebasecalled.fastq.tar.gz

echo "=== Extrayendo HG00733 ==="
tar -xzf HG00733.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > HG00733.fastq.gz
rm -rf projects/
rm -f HG00733.ONT.rebasecalled.fastq.tar.gz

echo "=== HG00733 listo ==="
ls -lh HG00733.fastq.gz
zcat HG00733.fastq.gz | head -4
du -sh ~/
