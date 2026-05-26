# RetroInspector-NF
Reimplementación en Nextflow DSL2 de [RetroInspector](https://github.com/javiercguard/retroinspector), un pipeline para 
la detección y caracterización de inserciones y deleciones de elementos transponibles (TEs) a partir de datos de 
secuenciación nanopore. Desarrollado como Trabajo de Fin de Máster en la Universidad de Murcia. Soporta los genomas de 
referencia **hg38** (GRCh38) y **T2T-CHM13v2.0**, con descarga automática de todos los recursos necesarios. ---
## Requisitos
- Nextflow >= 22.10 
- Conda o Mamba 
- Acceso a internet desde los nodos de cómputo (para descarga de recursos) 
- Más de 600 
GB de espacio en disco para una ejecución completa con T2T 

## Uso
Prepara un samplesheet CSV con tus muestras: 
sample_id,fastq 
HG00514,/ruta/HG00514.fastq.gz 
NA19240,/ruta/NA19240.fastq.gz 
...

Ejecuta con hg38 (descarga automática):

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir resultados \
    -profile dayhoff
```

Ejecuta con T2T-CHM13v2.0:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --outdir resultados \
    -profile dayhoff
```

Con tu propio genoma de referencia (sin anotación génica):

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /ruta/genoma.fa \
    --outdir resultados \
    -profile dayhoff
```

Para comparar pares de muestras, añade en `nextflow.config`:

```groovy
params {
    comparisons = [["HG00514", "NA19240"]]
}
```

---

## Parámetros

| Parámetro | Descripción | Por defecto |
|-----------|-------------|-------------|
| `--input` | Samplesheet CSV | null |
| `--reference` | Genoma de referencia propio | null |
| `--genome` | Genoma a descargar: `hg38` o `t2t` | hg38 |
| `--outdir` | Carpeta de resultados | results |
| `--mode` | `full` (con anotación) o `light` (sin anotación) | full |
| `--min_read_support` | Lecturas mínimas para aceptar una variante | 3 |
| `--threads` | Hilos de procesamiento | 32 |
| `--enrichment_pval` | P-valor para enriquecimiento funcional | 0.05 |
| `--dist_intra` | Distancia máxima para merge intra-muestra (bp) | 60 |
| `--dist_inter` | Distancia máxima para merge inter-muestra (bp) | 60 |
| `--comparisons` | Pares de muestras a comparar | [] |

---

## Salidas

results/
 ├── alns/ BAMs indexados por muestra
 ├── variants/ VCFs por caller, por muestra y conjuntos finales lax y strict
 ├── repeatmasker/ BED con familia de TE anotada, incluyendo SVA F1
 ├── rds/ Objetos R intermedios
 └── reports/ 
 ├── report.all.html Informe principal
 └── MUESTRA1_vs_MUESTRA2.html 

El informe HTML incluye estadísticas de filtrado, benchmark de genotyping,
distribución genómica de inserciones y deleciones, clasificación por familia
de TE, Manhattan plots de frecuencias alélicas, hallmarks de retrotransposición
(cola polyA y motivo de endonucleasa L1) y análisis de enriquecimiento funcional
(Gene Ontology, Disease Ontology y Network of Cancer Genes).

---

## Diferencias respecto al pipeline original

Esta reimplementación introduce las siguientes diferencias respecto al original.

El pipeline original está implementado en Snakemake. Esta versión usa Nextflow
DSL2 con módulos independientes, lo que facilita la ejecución en entornos HPC
con SLURM y permite el uso de `-resume` para reanudar ejecuciones interrumpidas.

**Soporte para T2T-CHM13v2.0**
El pipeline original solo soporta T2T en modo light, sin anotación génica. Esta
versión implementa el soporte completo en modo full, incorporando los siguientes
recursos específicos de T2T:

- Genoma CHM13v2.0 descargado desde el repositorio público del T2T Consortium
- BED de repeticiones chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed del mismo
  consorcio, con ordenación y reindexado previos con tabix --csi
- GTF de anotación génica hs1.ncbiRefSeq.gtf descargado desde UCSC
- Librería RepeatMasker construida con RepeatMasker 4.1.9 combinando la base
  de Dfam con el fichero humanAutoXYape.embl de Jessica Storer
  (github.com/jessicaStorer88/RepeatMasker_library_CHM13) más la partición 7
  de Dfam 3.9 (Mammalia), necesaria para detectar LINE, LTR y SVA
- Anotación génica construida desde el GTF con makeTxDbFromGFF de
  GenomicFeatures, dado que annotatr no tiene soporte nativo para CHM13

**Módulo HALLMARKS integrado**
En el pipeline original los hallmarks de retrotransposición son un análisis
independiente. Aquí están integrados como proceso Nextflow automático entre
R_GENOTYPING y R_REPORT, usando pysam y RapidFuzz.

**Migración de scripts**
Todos los scripts R han sido migrados de la API Snakemake (snakemake@input,
snakemake@output, snakemake@params) a commandArgs y parámetros de rmarkdown.
Los scripts Python han sido migrados de la API Snakemake a sys.argv.

**Otros cambios**
- Detección automática del formato BED CHM13 vs hg38 en getMeDeletions.py
- Manejo robusto de VCFs vacíos en assembly_alleles y merge
- Perfil de ejecución SLURM para el servidor dayhoff (cola eck-q)

---

## Cita

El trabajo en el que se basa:

Cuenca-Guardiola J, de la Morena-Barrio B, Corral J, Fernández-Breis JT.
Advanced analysis of retrotransposon variation in the human genome with
nanopore sequencing using RetroInspector. Scientific Reports. 2025;15:14489.
https://doi.org/10.1038/s41598-025-98847-7

---


Mónica Zapata Hernández
Máster Universitario en Bioinformática, UMU
