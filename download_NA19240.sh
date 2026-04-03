#!/bin/bash
#SBATCH --job-name=dl_NA19240
#SBATCH --output=logs/dl_NA19240_%j.log
#SBATCH --error=logs/dl_NA19240_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

cd ~/tfm/simulacion/input_data

echo "=== Descargando NA19240 ==="
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/NA19240.ONT.rebasecalled.fastq.tar.gz

echo "=== Extrayendo NA19240 ==="
tar -xzf NA19240.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > NA19240.fastq.gz
rm -rf projects/
rm -f NA19240.ONT.rebasecalled.fastq.tar.gz

echo "=== NA19240 listo ==="
ls -lh NA19240.fastq.gz
zcat NA19240.fastq.gz | head -4
du -sh ~/
