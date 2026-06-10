# RetroInspector-NF

Nextflow DSL2 reimplementation of [RetroInspector](https://github.com/javiercguard/retroinspector), a pipeline for the detection and characterization of transposable element (TE) insertions and deletions from Nanopore sequencing data.

This project was developed as a Master's Thesis at the University of Murcia. The pipeline keeps the original RetroInspector logic as closely as possible while extending the workflow to support both **hg38 / GRCh38** and **T2T-CHM13v2.0** reference genomes, and bundling all tooling into a single versioned Docker image.

---

## Requirements

* Nextflow (≥ 24.04)
* One of the following execution backends:
  * **Docker** (≥ 20.10) — recommended for reproducibility.
  * **Conda or Mamba** — fallback when Docker is not available.
* Internet access from compute nodes, if reference resources need to be downloaded automatically.
* Sufficient disk space for long-read alignments and intermediate files. Full T2T runs can require several hundred GB of storage depending on the number of samples.

---

## Input

Prepare a samplesheet CSV with one sample per row:

```csv
sample_id,fastq
HG00514,/path/HG00514.fastq.gz
NA19240,/path/NA19240.fastq.gz
HG00733,/path/HG00733.fastq.gz
```

---

## Usage

The pipeline supports two execution modes, selected via Nextflow profiles. The Docker profile is recommended; the Conda profile is kept as a fallback for environments where Docker is not available.

### Docker (recommended)

The Docker image `monicazh/retroinspector-nf:0.1.0` ships all four isolated conda environments (`retro-base`, `retro-sniffles2`, `retro-rm-t2t`, `retro-r`) needed by the pipeline. See the [Container architecture](#container-architecture) section for details.

```bash
# hg38 with automatic reference download
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir results_docker \
    -profile dayhoff,docker

# T2T-CHM13v2.0 with automatic reference download
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --outdir results_docker \
    -profile dayhoff,docker
```

### Conda fallback

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir results \
    -profile dayhoff,conda
```

### User-provided reference FASTA

For either execution mode, a custom reference can be supplied:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --reference /path/chm13v2.0.fa \
    --outdir results_docker \
    -profile dayhoff,docker
```

When `--reference` is provided, the pipeline uses that FASTA for alignment. The `--genome` parameter still determines which auxiliary resources are used, such as RepeatMasker reference annotations, gene annotations, and the T2T-specific RepeatMasker library.

---

## Container architecture

RetroInspector-NF is distributed as a single Docker image that bundles four isolated conda environments. This design avoids maintaining four separate images while keeping tool versions locked together at the image level.

### Image structure

| Conda environment    | Purpose                                                                        |
| -------------------- | ------------------------------------------------------------------------------ |
| `retro-base`         | Minimap2, samtools, cuteSV, mosdepth, surpyvor, bcftools, RepeatMasker (hg38). |
| `retro-sniffles2`    | Sniffles 2.x, isolated to avoid version conflicts with the legacy Sniffles.    |
| `retro-rm-t2t`       | RepeatMasker 4.1.9 + rmblast 2.14.1 for T2T-CHM13v2.0 annotation.              |
| `retro-r`            | R 4.2 + Bioconductor 3.16 (annotatr, clusterProfiler, rmarkdown, flextable).   |

### Building the image locally

```bash
docker build -t monicazh/retroinspector-nf:0.1.0 .
```

### Pulling the prebuilt image

```bash
docker pull monicazh/retroinspector-nf:0.1.0
```

### Container-level configuration

The Docker profile in `nextflow.config` includes:

* `-u $(id -u):$(id -g)` — runs the container as the host user, so files produced in the work directory remain owned by the user, not by root.
* `-e HOME=/tmp` — ensures non-root execution can write the micromamba shell caches required by the environments.

### Environment activation inside modules

Each process activates its target environment explicitly in its `script:` block by exporting `PATH` and invoking the binary by absolute path, for example:

```bash
R_ENV="/opt/conda/envs/retro-r"
export PATH="${R_ENV}/bin:${PATH}"
${R_ENV}/bin/Rscript ${projectDir}/bin/analysisPreparatory.R ...
```

This pattern is used because Nextflow overrides the container `ENTRYPOINT` with `/bin/bash` when launching processes, which bypasses any `activate-and-run.sh` wrapper embedded in the image. Activating the environment explicitly per process makes the pipeline independent of the container launch mechanism and reproducible on any Docker setup, with or without entrypoint customization.

---

## Pairwise comparisons

To compare pairs of samples, add the comparisons to `nextflow.config`:

```groovy
params {
    comparisons = [["HG00514", "NA19240"]]
}
```

---

## Parameters

| Parameter            | Description                                                | Default            |
| -------------------- | ---------------------------------------------------------- | ------------------ |
| `--input`            | Samplesheet CSV file                                       | `null`             |
| `--fastq_dir`        | Directory containing FASTQ/FASTQ.GZ files                  | `null`             |
| `--reference`        | Optional user-provided reference FASTA                     | `null`             |
| `--genome`           | Reference genome mode: `hg38` or `t2t`                     | `hg38`             |
| `--outdir`           | Output directory                                           | `results`          |
| `--mode`             | Analysis mode: `full` or `light`                           | `full`             |
| `--callers`          | Structural variant callers to use                          | `cuteSV,sniffles2` |
| `--min_read_support` | Minimum read support for variant calling                   | `3`                |
| `--threads`          | Number of processing threads                               | `32`               |
| `--enrichment_pval`  | P-value threshold for functional enrichment                | `0.05`             |
| `--dist_intra`       | Maximum distance for intra-sample insertion merging, in bp | `60`               |
| `--dist_inter`       | Maximum distance for inter-sample insertion merging, in bp | `60`               |

---

## Outputs

The output directory contains:

```text
results/
├── alns/              Indexed BAM files per sample
├── variants/          Caller-level and final VCF files
├── repeatmasker/      RepeatMasker-based TE annotation of insertions
├── rds/               Intermediate R objects
└── reports/           HTML reports
```

Main report files include:

```text
report.all.html
SAMPLE1_vs_SAMPLE2.html
```

The HTML report includes filtering statistics, genotyping summaries, genomic distribution of insertions and deletions, TE family classification, allele-frequency plots, retrotransposition hallmarks, and functional enrichment analyses.

---

## Relationship with the original RetroInspector pipeline

The original RetroInspector pipeline is implemented in Snakemake. RetroInspector-NF migrates the workflow to Nextflow DSL2 using modular processes, which improves execution on HPC systems with SLURM and enables interrupted runs to be resumed with `-resume`.

The main objective of this reimplementation is to preserve the original analytical logic while adapting the execution engine from Snakemake to Nextflow. Methodological changes are kept minimal unless required for the T2T extension or for container compatibility.

---

## hg38 / GRCh38 reference resources

For hg38, RetroInspector-NF follows the original RetroInspector reference logic.

The RepeatMasker reference BED is generated from the UCSC `hg38.fa.out.gz` file by:

1. downloading the RepeatMasker output from UCSC,
2. skipping the RepeatMasker header,
3. filtering TE classes matching `SINE`, `LINE`, `Retroposon`, `LTR`, and `DNA`,
4. extracting the same six columns used by the original pipeline,
5. compressing the BED with `bgzip`,
6. indexing it with `tabix --csi`.

This reproduces the original `rules/reference.smk` logic as closely as possible.

A non-overlapping TE reference file can also be generated with `bedtools merge`, following the original `prepare_te_nooverlap` rule. In the original pipeline, this file was generated as an auxiliary resource, but the active deletion annotation step used the indexed `repeatsReferenceTE.bed.gz` file. In RetroInspector-NF, any use of the non-overlapping version should therefore be treated as an optional extension rather than as the default original behavior.

---

## T2T-CHM13v2.0 extension

RetroInspector-NF extends the workflow to T2T-CHM13v2.0.

The T2T mode uses:

* the CHM13v2.0 reference genome,
* the CHM13 RepeatMasker BED annotation file `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed`,
* the `hs1.ncbiRefSeq.gtf` gene annotation from UCSC,
* a T2T-specific RepeatMasker library generated from the local RepeatMasker library and the `humanAutoXYape.embl` file from the CHM13 RepeatMasker library repository,
* Dfam 3.9 partitions 0 (root taxonomy) and 7 (Mammalia) for FamDB-based species lookup.

The CHM13 RepeatMasker BED is used as a reference annotation resource for deletion annotation. The custom RepeatMasker library is used separately during insertion annotation, when RepeatMasker is run on inserted sequences.

This distinction is important:

```text
GET_REFERENCE_REPEATS_T2T
    -> generates the CHM13 RepeatMasker BED used as a reference annotation resource

GET_REPEATMASKER_LIB_T2T
    -> prepares the T2T-specific RepeatMasker library used to annotate inserted sequences

GET_DFAM_PARTITION7
    -> downloads Dfam 3.9 partitions 0 and 7 required by famdb.py for human TE annotation
```

### Note on Dfam version

The pipeline pins the Dfam download to **release 3.9**, not to the `releases/current` alias. Dfam 4.0 (released May 2026) introduced FamDB schema 3.0.0, which is incompatible with the `famdb.py` version 2.0.0 bundled in RepeatMasker 4.1.9 used by this pipeline. Anchoring the download to a fixed release also protects the pipeline from future Dfam updates that would otherwise silently break the container by promoting incompatible files through `current`.

For Homo sapiens, two partitions are sufficient: partition 0 holds the taxonomy backbone (required by famdb.py for any species query), and partition 7 holds Mammalia.

### Note on RepeatMasker libdir construction

In T2T mode the REPEATMASKER process builds a private library directory (`rm_libdir/`) inside the work directory rather than writing symlinks into the read-only path `/opt/conda/envs/retro-rm-t2t/share/RepeatMasker/Libraries/famdb/` inside the image. The latter is not writable when the container runs as a non-root UID (the standard configuration in this pipeline), and the upstream-staged library also contains a non-directory entry at `famdb`. The private libdir is populated with symlinks to the upstream library entries plus the Dfam partition files, and is passed to RepeatMasker via its `-libdir` argument.

---

## Script migration

Several scripts from the original Snakemake workflow were adapted for execution inside Nextflow processes.

Python scripts were migrated from Snakemake-specific inputs and outputs to command-line arguments using `sys.argv` or `argparse`.

R scripts were migrated from the Snakemake API to `commandArgs` and explicit RMarkdown parameters where needed.

Examples of preserved or adapted scripts include:

* `getGoodAlts.py`
* `merge.py`
* `genotype.py`
* `vcfToBedForMosdepth.py`
* `checkInsertionsMultiSample.py`
* `getMeDeletions.py`
* `svaf.py`
* `analysisPreparatory.R`
* `analysisGenotyping.R`
* `enrichment.R`
* `report.Rmd`
* `comparison.Rmd`

A new `hallmarks.py` script was added for the retrotransposition hallmark detection module.

---

## Main differences from the original pipeline

RetroInspector-NF introduces the following changes relative to the original RetroInspector workflow:

### Workflow engine

* Original pipeline: Snakemake.
* Current implementation: Nextflow DSL2.

### Execution model

* Modular Nextflow processes.
* SLURM-compatible execution profile.
* Support for `-resume`.

### Containerization

* Single versioned Docker image (`monicazh/retroinspector-nf:0.1.0`) bundling four isolated conda environments.
* Docker run options configured for non-root user execution.

### Reference support

* Original pipeline: hg38 / GRCh38.
* RetroInspector-NF: hg38 / GRCh38 and T2T-CHM13v2.0.

### T2T-specific additions

* CHM13v2.0 reference genome support.
* CHM13 RepeatMasker BED support.
* hs1 gene annotation support.
* T2T-specific RepeatMasker library preparation.
* Adapted parsing of hg38 and CHM13 RepeatMasker BED formats where required.
* T2T-specific VCF header derived from CHM13v2.0 contig lengths.

### Data layer pinning

* Dfam download anchored to release 3.9 to guarantee FamDB schema compatibility with the bundled RepeatMasker version, regardless of upstream alias changes.

### Added analysis modules

* Optional retrotransposition hallmark detection, including polyA/polyT tails and L1 endonuclease motif detection.
* T2T-specific reporting components.

---

## Citation

This work is based on the original RetroInspector pipeline:

Cuenca-Guardiola J, de la Morena-Barrio B, Corral J, Fernández-Breis JT.
Advanced analysis of retrotransposon variation in the human genome with nanopore sequencing using RetroInspector.
*Scientific Reports*. 2025;15:14489.
https://doi.org/10.1038/s41598-025-98847-7

---

Mónica Zapata Hernández — University of Murcia
