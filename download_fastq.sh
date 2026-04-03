#!/bin/bash
#SBATCH --job-name=download_hgsvc
#SBATCH --output=logs/download_%j.log
#SBATCH --error=logs/download_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

cd ~/tfm/simulacion/input_data

# Descargar los tres tar.gz
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00514.ONT.rebasecalled.fastq.tar.gz
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/HG00733.ONT.rebasecalled.fastq.tar.gz
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/NA19240.ONT.rebasecalled.fastq.tar.gz

# Extraer y concatenar HG00514
tar -xzf HG00514.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > HG00514.fastq.gz
rm -rf projects/
echo "HG00514 listo: $(ls -lh HG00514.fastq.gz)"

# Extraer y concatenar HG00733
tar -xzf HG00733.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > HG00733.fastq.gz
rm -rf projects/
echo "HG00733 listo: $(ls -lh HG00733.fastq.gz)"

# Extraer y concatenar NA19240
tar -xzf NA19240.ONT.rebasecalled.fastq.tar.gz
find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > NA19240.fastq.gz
rm -rf projects/
echo "NA19240 listo: $(ls -lh NA19240.fastq.gz)"

echo "Todo completado"
ls -lh
