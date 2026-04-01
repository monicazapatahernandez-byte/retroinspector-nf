# RetroInspector-NF

Pipeline en Nextflow para la detección y caracterización de inserciones y 
deleciones de retrotransposones a partir de datos de secuenciación por nanoporos.

Este pipeline es una reimplementación en Nextflow del pipeline RetroInspector 
(originalmente en Snakemake), desarrollado como Trabajo de Fin de Máster en la 
Universidad de Murcia.

## Pipeline original
https://github.com/javiercguard/retroinspector

## Uso

### Con genoma hg38 propio
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /ruta/hg38.fa \
    --outdir resultados
```

### Descarga automática de hg38
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir resultados
```

### Con genoma T2T (modo light)
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --outdir resultados
```

## Parámetros principales

| Parámetro | Descripción | Por defecto |
|-----------|-------------|-------------|
| `--input` | Samplesheet CSV con columnas sample_id y fastq | null |
| `--reference` | Ruta a genoma de referencia propio | null |
| `--genome` | Genoma a usar: hg38 o t2t | hg38 |
| `--outdir` | Directorio de resultados | results |
| `--mode` | Modo de ejecución: full o light | full |
| `--min_read_support` | Mínimo de lecturas de soporte | 3 |
| `--threads` | Número de hilos | 32 |

## Samplesheet

Archivo CSV con este formato:
```
sample_id,fastq
HG00514,/ruta/HG00514.fastq.gz
HG00733,/ruta/HG00733.fastq.gz
NA19240,/ruta/NA19240.fastq.gz
```

## Autor
Mónica Zapata Hernández  

