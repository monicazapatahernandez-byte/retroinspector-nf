# RetroInspector-NF

Reimplementacion en Nextflow del pipeline RetroInspector, desarrollada como
Trabajo de Fin de Master en la Universidad de Murcia. Detecta y caracteriza
inserciones y deleciones de elementos transponibles (TEs) a partir de datos
de secuenciacion por nanoporos.

Pipeline original (Snakemake): https://github.com/javiercguard/retroinspector

---

## Que hace este pipeline?

Los elementos transponibles son secuencias de ADN que pueden moverse por el
genoma. Suponen el 45% del genoma humano y tienen implicaciones en
enfermedades, regulacion genica y variabilidad genetica. Detectarlos con
nanoporos es mas preciso que con lecturas cortas.

El pipeline recibe FASTQs de nanopore y un genoma de referencia, y produce
un informe HTML con las inserciones y deleciones de TEs, su anotacion genica
y analisis de enriquecimiento.

---

## Pasos del pipeline

  1. Alineamiento        Minimap2 contra el genoma de referencia
  2. Variant calling     CuteSV y Sniffles2 en paralelo
  3. Ensamblaje          Reconstruccion de la secuencia insertada
  4. Merge intra-muestra Fusion de llamadas de ambos callers por muestra
  5. Genotyping          Cobertura con Mosdepth y calculo de genotipos
  6. Merge inter-muestra Fusion de todas las muestras
  7. RepeatMasker        Identificacion de TEs y deteccion de SVA F1 con BLAST
  8. Deleciones          Deteccion de deleciones de TE frente al genoma
  9. Analisis R          Anotacion, enriquecimiento e informe HTML

El paso 9 solo se ejecuta en modo full (por defecto). En modo light se omite
la anotacion, permitiendo usar genomas alternativos como T2T.

---

## Uso

Formato del samplesheet (CSV):

  sample_id,fastq
  HG00514,/ruta/HG00514.fastq.gz
  HG00733,/ruta/HG00733.fastq.gz
  NA19240,/ruta/NA19240.fastq.gz

Con genoma propio:

  nextflow run main.nf \
      --input samplesheet.csv \
      --reference /ruta/hg38.fa \
      --outdir resultados \
      -profile dayhoff

Con descarga automatica de hg38:

  nextflow run main.nf \
      --input samplesheet.csv \
      --genome hg38 \
      --outdir resultados \
      -profile dayhoff

Con genoma T2T (modo light, sin anotacion):

  nextflow run main.nf \
      --input samplesheet.csv \
      --genome t2t \
      --outdir resultados \
      -profile dayhoff

---

## Parametros principales

  --input              Samplesheet CSV                         (null)
  --reference          Genoma de referencia propio             (null)
  --genome             Genoma a descargar: hg38 o t2t          (hg38)
  --outdir             Carpeta de resultados                   (results)
  --mode               full o light                            (full)
  --min_read_support   Lecturas minimas para aceptar variante  (3)
  --threads            Hilos de procesamiento                  (32)
  --all_prefix         Prefijo de ficheros de salida           (all)
  --enrichment_pval    P-valor para enriquecimiento            (0.05)
  --dist_intra         Distancia merge intra-muestra (bp)      (60)
  --dist_inter         Distancia merge inter-muestra (bp)      (60)
  --comparisons        Pares de muestras a comparar            ([])

Para comparar muestras, añadir en nextflow.config:

  params {
      comparisons = [["HG00514", "NA19240"]]
  }

---

## Resultados

  results/
  |-- alns/          Alineamientos BAM
  |-- variants/      VCFs por muestra y conjuntos finales
  |-- repeatmasker/  Anotacion RepeatMasker con SVA F1
  |-- rds/           Objetos R intermedios
  |-- reports/
      |-- report.all.html            Informe principal
      |-- MUESTRA1_vs_MUESTRA2.html  Informe de comparacion

El informe HTML incluye:
  - Estadisticas de filtrado (criterio lax y strict)
  - Benchmark de genotyping (falsos positivos y negativos)
  - Distribucion genomica de inserciones y deleciones
  - Clasificacion por familia de TE (SVA, Alu, L1 detallados)
  - Manhattan plots de frecuencias alelicas
  - Analisis de enriquecimiento (GO, Disease Ontology, NCG)

---

## Diferencias respecto al pipeline original

  - Reimplementado completamente en Nextflow DSL2
  - Scripts R migrados de API Snakemake a parametros de rmarkdown
  - Deteccion de SVA F1 integrada en el proceso RepeatMasker
  - Soporte para comparacion entre pares de muestras (R_COMPARE)
  - Soporte para genoma T2T en modo light

---

## Autor

Monica Zapata Hernandez

