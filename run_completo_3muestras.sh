#!/bin/bash
#SBATCH --job-name=retro_3muestras
#SBATCH --output=logs/retro_3muestras_%j.log
#SBATCH --error=logs/retro_3muestras_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/alumno27/miniconda3/envs/nf-core

DATA_DIR=~/tfm/retroinspector-nf/input_data
PIPE_DIR=~/tfm/retroinspector-nf
BASE_URL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled

echo "=== INICIO: $(date) ==="
du -sh ~

# ============================================================
# Función para procesar una muestra
# ============================================================
procesar_muestra() {
    local SAMPLE=$1
    local TARFILE="${SAMPLE}.ONT.rebasecalled.fastq.tar.gz"
    local FASTQFILE="${SAMPLE}.fastq.gz"

    echo ""
    echo "=== [$SAMPLE] INICIO: $(date) ==="
    du -sh ~

    if [ -f "$DATA_DIR/$FASTQFILE" ]; then
        echo "=== [$SAMPLE] $FASTQFILE ya existe, saltando descarga ==="
    else
        rm -f $DATA_DIR/$TARFILE
        echo "=== [$SAMPLE] Descargando ==="
        cd $DATA_DIR
        wget -c ${BASE_URL}/${TARFILE}

        echo "=== [$SAMPLE] Verificando integridad ==="
        if ! gzip -t $TARFILE 2>/dev/null; then
            echo "ERROR [$SAMPLE]: tar.gz corrupto"
            exit 1
        fi

        echo "=== [$SAMPLE] Extrayendo ==="
        tar -xzf $TARFILE
        find projects/ -name "*.fastq.gz" | sort | xargs zcat | gzip > $FASTQFILE
        rm -rf projects/
        rm -f $TARFILE

        echo "=== [$SAMPLE] Verificando FASTQ ==="
        N_READS=$(zcat $FASTQFILE | head -400 | grep "^@" | wc -l)
        echo "Reads en primeras 400 líneas: $N_READS"
        if [ "$N_READS" -eq 0 ]; then
            echo "ERROR [$SAMPLE]: FASTQ vacío o corrupto, abortando"
            exit 1
        fi
    fi

    echo "Tamaño FASTQ: $(du -sh $DATA_DIR/$FASTQFILE | cut -f1)"
    echo "=== [$SAMPLE] FASTQ listo ==="
    du -sh ~
}

# ============================================================
# Descargar y preparar las 3 muestras
# ============================================================
procesar_muestra HG00514
procesar_muestra HG00733
procesar_muestra NA19240

# ============================================================
# Lanzar Nextflow con las 3 muestras
# ============================================================
echo ""
echo "=== Lanzando Nextflow con 3 muestras: $(date) ==="
cd $PIPE_DIR

/home/alumno27/miniconda3/envs/nf-core/bin/nextflow run main.nf \
    --input test_data/samplesheet.csv \
    --reference /home/software/nanopore/hg38.fa \
    --outdir results_3muestras \
    --threads 32 \
    -profile dayhoff \
    -work-dir work_3muestras \
    -resume

echo "=== FIN: $(date) ==="
du -sh ~
