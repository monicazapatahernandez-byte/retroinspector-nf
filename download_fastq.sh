#!/bin/bash
#SBATCH --job-name=download_hgsvc
#SBATCH --output=logs/download_%j.log
#SBATCH --error=logs/download_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=eck-q

cd ~/tfm/simulacion/input_data

# Descargar solo HG00514 primero para validar el pipeline
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00514.ONT.rebasecalled.fastq.tar.gz

# Extraer según instrucciones del repositorio
tar -zxf HG00514.ONT.rebasecalled.fastq.tar.gz --wildcards '*.fastq*' -O 2>/dev/null | gzip > HG00514.fastq.gz

echo "Descarga y extracción completadas"
ls -lh ~/tfm/simulacion/input_data/
