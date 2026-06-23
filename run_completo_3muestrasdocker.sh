#!/bin/bash
# RetroInspector-NF — Other servers
#
# Normal:
#   ./run_completo_3muestrasdocker.sh
# Closed SSH:
#   nohup ./run_completo_3muestrasdocker.sh > logs/retro_hg38.log 2> logs/retro_hg38.err &


PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PIPE_DIR}/input_data"
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled"

mkdir -p "${DATA_DIR}" "${PIPE_DIR}/logs"

echo "=== INICIO: $(date) ==="
nextflow -version
docker --version

procesar_muestra() {
    local SAMPLE=$1
    local TARFILE="${SAMPLE}.ONT.rebasecalled.fastq.tar.gz"
    local FASTQFILE="${SAMPLE}.fastq.gz"

    if [ -f "${DATA_DIR}/${FASTQFILE}" ]; then
        echo "[$SAMPLE] FASTQ ya existe, saltando"
        return 0
    fi

    echo "[$SAMPLE] Descargando..."
    wget -c "${BASE_URL}/${TARFILE}" -O "${DATA_DIR}/${TARFILE}"

    echo "[$SAMPLE] Extrayendo..."
    tar -xzf "${DATA_DIR}/${TARFILE}" -C "${DATA_DIR}"

    echo "[$SAMPLE] Consolidando FASTQs..."
    find "${DATA_DIR}" -name "*.fastq.gz" ! -name "${FASTQFILE}" | sort \
        | xargs zcat | gzip > "${DATA_DIR}/${FASTQFILE}"

    find "${DATA_DIR}" -mindepth 1 -maxdepth 1 -type d | xargs rm -rf
    rm -f "${DATA_DIR}/${TARFILE}"
    echo "[$SAMPLE] Listo"
}

procesar_muestra HG00514
procesar_muestra NA19240
procesar_muestra HG00733

cd "${PIPE_DIR}"
nextflow run main.nf \
    --input test_data/samplesheet.csv \
    --genome hg38 \
    --mode full \
    --outdir results_dockerhg38 \
    -profile docker \
    -work-dir work_dockerhg38   \
    -with-dag reports/dag_hg38.html \
    -resume


echo "=== FIN: $(date) ==="
