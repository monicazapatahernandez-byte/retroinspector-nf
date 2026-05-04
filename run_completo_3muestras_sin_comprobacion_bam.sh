#!/bin/bash
#SBATCH --job-name=retro_3muestras
#SBATCH --output=logs/retro_3muestras_%j.log
#SBATCH --error=logs/retro_3muestras_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=eck-q

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/alumno27/miniconda3/envs/nf-core

DATA_DIR=/home/alumno27/tfm/retroinspector-nf/input_data
PIPE_DIR=/home/alumno27/tfm/retroinspector-nf
BASE_URL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled

echo "=== INICIO: $(date) ==="
du -sh ~

procesar_muestra() {
    local SAMPLE=$1
    local TARFILE="${SAMPLE}.ONT.rebasecalled.fastq.tar.gz"
    local FASTQFILE="${SAMPLE}.fastq.gz"

    echo ""
    echo "=== [$SAMPLE] INICIO: $(date) ==="
    du -sh ~

    if [ -f "${DATA_DIR}/${FASTQFILE}" ]; then
        echo "=== [$SAMPLE] $FASTQFILE ya existe, saltando descarga ==="
    else
        rm -f "${DATA_DIR}/${TARFILE}"
        echo "=== [$SAMPLE] Descargando ==="
        wget -c "${BASE_URL}/${TARFILE}" -O "${DATA_DIR}/${TARFILE}"

        echo "=== [$SAMPLE] Verificando integridad del tar ==="
        if ! tar -tzf "${DATA_DIR}/${TARFILE}" > /dev/null 2>&1; then
            echo "ERROR [$SAMPLE]: tar.gz corrupto"
            exit 1
        fi

        echo "=== [$SAMPLE] Extrayendo ==="
        tar -xzf "${DATA_DIR}/${TARFILE}" -C "${DATA_DIR}"

        echo "=== [$SAMPLE] Buscando ficheros fastq ==="
        FASTQ_FILES=$(find "${DATA_DIR}" -name "*.fastq.gz" ! -name "${FASTQFILE}" | sort)
        if [ -z "$FASTQ_FILES" ]; then
            FASTQ_FILES=$(find "${DATA_DIR}" -name "*.fastq" | sort)
            if [ -z "$FASTQ_FILES" ]; then
                echo "ERROR [$SAMPLE]: no se encontraron ficheros fastq tras extracción"
                exit 1
            fi
            echo "Encontrados $(echo "$FASTQ_FILES" | wc -l) ficheros .fastq"
            echo "$FASTQ_FILES" | xargs cat | gzip > "${DATA_DIR}/${FASTQFILE}"
        else
            echo "Encontrados $(echo "$FASTQ_FILES" | wc -l) ficheros .fastq.gz"
            echo "$FASTQ_FILES" | xargs zcat | gzip > "${DATA_DIR}/${FASTQFILE}"
        fi

        # Limpiar carpetas extraídas
        find "${DATA_DIR}" -mindepth 1 -maxdepth 1 -type d | xargs rm -rf
        rm -f "${DATA_DIR}/${TARFILE}"

        echo "=== [$SAMPLE] Verificando FASTQ ==="
        if ! gzip -t "${DATA_DIR}/${FASTQFILE}" 2>&1; then
            echo "ERROR [$SAMPLE]: FASTQ corrupto tras extracción"
            exit 1
        fi
        N_READS=$(zcat "${DATA_DIR}/${FASTQFILE}" | head -4000 | awk 'NR%4==1 && /^@/{count++} END{print count}')
        echo "Reads en primeras 4000 líneas: $N_READS"
        if [ "$N_READS" -lt 10 ]; then
            echo "ERROR [$SAMPLE]: FASTQ sospechoso (menos de 10 reads en 4000 líneas)"
            exit 1
        fi
    fi

    echo "Tamaño FASTQ: $(du -sh ${DATA_DIR}/${FASTQFILE} | cut -f1)"
    echo "=== [$SAMPLE] FASTQ listo ==="
    du -sh ~
}

procesar_muestra HG00514
procesar_muestra HG00733
procesar_muestra NA19240

echo ""
echo "=== Lanzando Nextflow con 3 muestras: $(date) ==="
cd "$PIPE_DIR"
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
