# RetroInspector-NF

Nextflow DSL2 reimplementation of RetroInspector, a pipeline for detection and characterization of transposable element (TE) insertions and deletions from Nanopore sequencing data.

This project was developed as a Master's Thesis (TFM) at the University of Murcia (UMU). The pipeline keeps the original RetroInspector logic as closely as possible while extending the workflow to support both hg38 / GRCh38 and T2T-CHM13v2.0 reference genomes, and bundling all tooling into a single versioned Docker image.

## 1. Installation

RetroInspector-NF runs directly from a clone of this repository. The host system only needs:

- Nextflow ≥ 24.04 (validated on 25.10.4). The pipeline is DSL2 and uses features only available on recent Nextflow releases.
- Docker ≥ 20.10. Docker is the recommended (and primary) execution backend. The pipeline ships a single image that bundles all required tools, so nothing else needs to be installed on the host.
- Java ≥ 17 on the host that launches Nextflow.
- Conda / Mamba (optional). Kept as a fallback for environments where Docker is not available; see Section 3.3.

The compute nodes that execute the pipeline need outbound Internet access to:

- pull the Docker image `monicazh/retroinspector-nf:0.1.0` from Docker Hub on the first run,
- download reference resources (FASTA, RepeatMasker BED, gene annotation, Dfam partitions) the first time a given `--genome` mode is launched, unless the user provides them locally via `--reference`.

Once the image and the references have been cached, the pipeline can be re-executed offline within the same work directory.

For resource planning, a full run with three HGSV samples on both references needs around 1 TB of free disk and at least 64 GB of RAM on the execution node (MINIMAP2_ALIGN is configured to request 48 GB and REPEATMASKER 32 GB; both can be lowered in `nextflow.config` for smaller machines). Disk usage breaks down roughly as:

- ~250–350 GB for raw and consolidated FASTQ inputs,
- ~300–500 GB for the Nextflow work directory (alignments and intermediates),
- ~50–100 GB for cached reference resources (hg38 + T2T-CHM13v2.0, RepeatMasker libraries, Dfam partitions),
- the rest for the published outputs under `results_docker*/`.

To get the code:

```bash
git clone https://github.com/monicazapatahernandez-byte/retroinspector-nf.git
cd retroinspector-nf
```

There is nothing else to install at this stage. The first call to `nextflow run` will pull the Docker image automatically.

## 2. Configuration

The pipeline is configured through two files at the root of the repository: `main.nf` (the workflow itself) and `nextflow.config` (parameters, profiles and per-process resources). In normal use only `nextflow.config` needs to be edited; `main.nf` is touched only when the workflow logic changes.

### 2.1 main.nf

`main.nf` imports the DSL2 modules under `modules/` and connects them into the full workflow: sample ingestion from a CSV samplesheet or a FASTQ directory, reference resolution (hg38 or T2T-CHM13v2.0, either user-provided via `--reference` or auto-downloaded), alignment with Minimap2, SV calling with CuteSV and Sniffles2, allele assembly, intrapatient and interpatient merging, insertion and deletion genotyping with Mosdepth, RepeatMasker annotation of inserted sequences, R-based downstream analysis (preparation, genotyping, enrichment, comparison, report), and hallmark detection (polyA tails and L1 endonuclease motif).

`main.nf` does not need to be edited to change samples, references or analysis options: everything is exposed through `params` in `nextflow.config` and through command-line flags.

### 2.2 nextflow.config

`nextflow.config` defines three blocks: pipeline parameters, execution profiles, and per-process resources.

**Pipeline parameters.** All adjustable parameters live under the `params { ... }` block and can be overridden from the command line. The defaults shipped with the repository are:

| Parameter            | Description                                              | Default            |
|----------------------|----------------------------------------------------------|--------------------|
| `--input`            | Samplesheet CSV file                                     | `null`             |
| `--fastq_dir`        | Directory containing FASTQ/FASTQ.GZ files                | `null`             |
| `--reference`        | Optional user-provided reference FASTA                   | `null`             |
| `--genome`           | Reference genome mode: `hg38` or `t2t`                   | `null`             |
| `--outdir`           | Output directory                                         | `results`          |
| `--mode`             | Analysis mode: `full` or `light`                         | `null`             |
| `--callers`          | SV callers to use (comma-separated, no spaces)           | `cuteSV,sniffles2` |
| `--min_read_support` | Minimum read support for variant calling                 | `3`                |
| `--threads`          | Number of processing threads                             | `32`               |
| `--enrichment_pval`  | P-value threshold for functional enrichment              | `0.05`             |
| `--dist_intra`       | Max distance for intra-sample insertion merging, in bp   | `60`               |
| `--dist_inter`       | Max distance for inter-sample insertion merging, in bp   | `60`               |

`--genome` and `--mode` must be set explicitly at run time; the included execution scripts already do so. `--mode full` runs the complete workflow including the R analyses and the HTML report; `--mode light` stops after VCF generation and skips all downstream R steps, which is useful for quick variant-calling-only runs.

**Pairwise comparisons.** Pairwise sample comparisons are configured directly in `nextflow.config`, since list-of-lists values are awkward to pass on the command line:

```groovy
params {
    comparisons = [["HG00514", "NA19240"]]
}
```

Each inner list is a pair of `sample_id` values present in the samplesheet. Several pairs can be added to the outer list.

**Execution profiles.** Two profiles are defined:

- `docker` — the recommended (and primary) profile. Enables Docker and sets the pipeline container to `monicazh/retroinspector-nf:0.1.0`. It also passes `-u $(id -u):$(id -g)` so files produced in the work directory remain owned by the host user.
- `conda` — fallback profile for environments without Docker. Enables Conda with Mamba and caches environments under `${projectDir}/.conda_envs`.

Only one profile is used per run, selected with the `-profile` flag.

**Per-process resources.** The `process { ... }` block sets CPU and memory budgets per module (Minimap2, CuteSV, Sniffles2, RepeatMasker, R steps, T2T RepeatMasker library preparation) and a retry strategy with up to two retries on transient failures. These values are tuned for an HPC-like node; they can be lowered for laptops or workstations by editing the block.

### 2.3 Container-level configuration

When the `docker` profile is active, each container runs with `-u $(id -u):$(id -g)` so files produced in the work directory remain owned by the host user (not by root).

Each Nextflow process activates its target conda environment explicitly inside its `script:` block by exporting `PATH` and invoking the binary by absolute path, for example:

```bash
R_ENV="/opt/conda/envs/retro-r"
export PATH="${R_ENV}/bin:${PATH}"
${R_ENV}/bin/Rscript ${projectDir}/bin/analysisPreparatory.R ...
```

This pattern is used because the Dockerfile does not define a custom `ENTRYPOINT`. Activating the environment per process makes the pipeline independent of container entrypoint behavior and reproducible on any Docker setup.

## 3. Execution

The pipeline is executed with `nextflow run main.nf` plus the appropriate flags and profile. The repository ships ready-to-use Bash scripts that automate the full execution on the HGSV samples (HG00514, NA19240, HG00733) for both reference modes, plus SLURM submission scripts under `scripts/` for HPC use. They can be used as-is on similar setups, or adapted to other clusters and other sample sets.

### 3.1 Samplesheet

Regardless of the execution mode, the samplesheet is a CSV with one sample per row:

```csv
sample_id,fastq
HG00514,/path/HG00514.fastq.gz
NA19240,/path/NA19240.fastq.gz
HG00733,/path/HG00733.fastq.gz
```

A reference samplesheet is included under `test_data/samplesheet.csv`. Alternatively, `--fastq_dir` can point the pipeline at a directory of FASTQ/FASTQ.GZ files; in that case the `sample_id` is taken from each file's `simpleName`.

### 3.2 Test data

The included scripts download the validation samples (HG00514, NA19240, HG00733) from the 1000 Genomes Project FTP, specifically the rebasecalled Nanopore tarballs at:

```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181210_ONT_rebasecalled/
```

Each sample is distributed as `<SAMPLE>.ONT.rebasecalled.fastq.tar.gz`. The scripts extract each archive, concatenate the per-flowcell FASTQs into a single `<SAMPLE>.fastq.gz` under `input_data/`, and verify integrity before launching the pipeline. To use other samples, just place their FASTQs under `input_data/` (or any other directory) and point the samplesheet at them.

### 3.3 Docker (recommended)

Docker is the primary execution mode. The image `monicazh/retroinspector-nf:0.1.0` is pulled automatically from Docker Hub on the first run.

Two wrapper scripts are provided at the root of the repository. Both first download and consolidate the HGSV FASTQ files, then launch the pipeline with `-profile docker`:

- `run_completo_3muestrasdocker.sh` — hg38 mode.
- `run_completo_3muestrast2tdocker.sh` — T2T-CHM13v2.0 mode.

They were used during development and validation on the dayhoff HPC at the University of Murcia, but they only depend on a standard Docker installation and a working Nextflow, so they are intended to be reusable on other systems as well.

Run them directly:

```bash
./run_completo_3muestrasdocker.sh        # hg38
./run_completo_3muestrast2tdocker.sh     # T2T-CHM13v2.0
```

Or in the background with SSH detached:

```bash
nohup ./run_completo_3muestrasdocker.sh    > logs/retro_hg38.log 2> logs/retro_hg38.err &
nohup ./run_completo_3muestrast2tdocker.sh > logs/retro_t2t.log  2> logs/retro_t2t.err  &
```

For SLURM clusters, equivalent submission scripts are provided under `scripts/` (for example `scripts/dayhoff_hg38.sh`). They wrap the same logic into an `sbatch`-ready file and write logs under `logs/`.

Each script writes its work directory and results into a separate folder (`work_dockerhg38` / `results_dockerhg38` for hg38, `work_dockert2t` / `results_dockert2t` for T2T) so the two modes can coexist on the same checkout without overwriting each other.

For other sample sets or other servers, the same execution can be launched directly with Nextflow, without the wrapping script:

```bash
# hg38 with automatic reference download
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --mode full \
    --outdir results_dockerhg38 \
    -profile docker \
    -work-dir work_dockerhg38 \
    -resume

# T2T-CHM13v2.0 with automatic reference download
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --mode full \
    --outdir results_dockert2t \
    -profile docker \
    -work-dir work_dockert2t \
    -resume
```

A user-provided reference FASTA can be passed in either mode:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome t2t \
    --reference /path/chm13v2.0.fa \
    --mode full \
    --outdir results_dockert2t \
    -profile docker
```

When `--reference` is provided, the pipeline uses that FASTA for alignment. The `--genome` parameter still determines which auxiliary resources are used: RepeatMasker reference annotations, gene annotations, and the T2T-specific RepeatMasker library.

### 3.4 Conda fallback

The `conda` profile is kept as a fallback for environments where Docker is not available. It builds the four required environments locally with Mamba and caches them under `${projectDir}/.conda_envs`. It is not the recommended execution mode and is not exercised by the included scripts:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --mode full \
    --outdir results \
    -profile conda
```

## 4. Container architecture

RetroInspector-NF is distributed as a single Docker image that bundles four isolated conda environments. This design avoids maintaining four separate images while keeping tool versions locked together at the image level.

| Conda environment   | Purpose                                                                          |
|---------------------|----------------------------------------------------------------------------------|
| `retro-base`        | Minimap2, samtools, cuteSV, mosdepth, surpyvor, bcftools, RepeatMasker (hg38).   |
| `retro-sniffles2`   | Sniffles 2.x, isolated to avoid version conflicts with legacy Sniffles.          |
| `retro-rm-t2t`      | RepeatMasker 4.1.9 + rmblast 2.14.1 for T2T-CHM13v2.0 annotation.                |
| `retro-r`           | R 4.2 + Bioconductor 3.16 (annotatr, clusterProfiler, rmarkdown, flextable).     |

The image can be rebuilt from the `Dockerfile` and the four `*.yaml` env files at the root of the repository, which fully pin the toolchain:

```bash
# Build locally from the repository (reproduces the exact image)
docker build -t monicazh/retroinspector-nf:0.1.0 .

# Or pull the prebuilt image from Docker Hub
docker pull monicazh/retroinspector-nf:0.1.0
```

## 5. Outputs

The output directory is named after `--outdir`. The included scripts use `results_dockerhg38/` for hg38 runs and `results_dockert2t/` for T2T runs. Inside, the layout is:

```
results_docker<mode>/
├── alns/              Indexed BAM files per sample
├── data/              Cached reference resources (FASTA, BED, GTF, Dfam, RM library)
├── variants/          Caller-level and final VCF files
├── repeatmasker/      RepeatMasker-based TE annotation of insertions
├── rds/               Intermediate R objects
├── reports/           HTML reports
└── pipeline_info/     Nextflow timeline and execution reports
```

Main report files include `report.all.html` and one `SAMPLE1_vs_SAMPLE2.html` per configured comparison. The HTML report contains filtering statistics, genotyping summaries, genomic distribution of insertions and deletions, TE family classification, allele-frequency plots, retrotransposition hallmarks, and functional enrichment analyses.

## 6. Relationship with the original RetroInspector pipeline

The original RetroInspector pipeline is implemented in Snakemake. RetroInspector-NF migrates the workflow to Nextflow DSL2 using modular processes, which improves execution on HPC systems with SLURM and enables interrupted runs to be resumed with `-resume`.

The main objective of this reimplementation is to preserve the original analytical logic while adapting the execution engine from Snakemake to Nextflow. Methodological changes are kept minimal unless required for the T2T extension or for container compatibility.

## 7. hg38 / GRCh38 reference resources

For hg38, RetroInspector-NF follows the original RetroInspector reference logic. The RepeatMasker reference BED is generated from the UCSC `hg38.fa.out.gz` file by downloading the RepeatMasker output from UCSC, skipping the RepeatMasker header, filtering TE classes matching SINE, LINE, Retroposon, LTR, and DNA, extracting the same six columns used by the original pipeline, compressing the BED with `bgzip`, and indexing it with `tabix --csi`.

A non-overlapping TE reference file can also be generated with `bedtools merge`, following the original `prepare_te_nooverlap` rule. In the original pipeline, this file was generated as an auxiliary resource, but the active deletion annotation step used the indexed `repeatsReferenceTE.bed.gz` file. In RetroInspector-NF, any use of the non-overlapping version should therefore be treated as an optional extension rather than as the default original behavior.

## 8. T2T-CHM13v2.0 extension

RetroInspector-NF extends the workflow to T2T-CHM13v2.0. The T2T mode uses:

- the CHM13v2.0 reference genome,
- the CHM13 RepeatMasker BED annotation file `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed`,
- the `hs1.ncbiRefSeq.gtf` gene annotation from UCSC,
- a T2T-specific RepeatMasker library generated from the local RepeatMasker library and the `humanAutoXYape.embl` file from the CHM13 RepeatMasker library repository,
- Dfam 3.9 partitions 0 (root taxonomy) and 7 (Mammalia) for FamDB-based species lookup.

The CHM13 RepeatMasker BED is used as a reference annotation resource for deletion annotation. The custom RepeatMasker library is used separately during insertion annotation, when RepeatMasker is run on inserted sequences:

```
GET_REFERENCE_REPEATS_T2T
    -> generates the CHM13 RepeatMasker BED used as a reference annotation resource

GET_REPEATMASKER_LIB_T2T
    -> prepares the T2T-specific RepeatMasker library used to annotate inserted sequences

GET_DFAM_PARTITION7
    -> downloads Dfam 3.9 partitions 0 and 7 required by famdb.py for human TE annotation
```

**Note on Dfam version.** The pipeline pins the Dfam download to release 3.9, not to the `releases/current` alias. Dfam 4.0 (released May 2026) introduced FamDB schema 3.0.0, which is incompatible with the `famdb.py` version 2.0.0 bundled in RepeatMasker 4.1.9 used by this pipeline. Anchoring the download to a fixed release also protects the pipeline from future Dfam updates that would otherwise silently break the container by promoting incompatible files through `current`. For *Homo sapiens*, two partitions are sufficient: partition 0 holds the taxonomy backbone (required by `famdb.py` for any species query), and partition 7 holds Mammalia.

**Note on RepeatMasker libdir construction.** In T2T mode the `REPEATMASKER` process builds a private library directory (`rm_libdir/`) inside the work directory rather than writing symlinks into the read-only path `/opt/conda/envs/retro-rm-t2t/share/RepeatMasker/Libraries/famdb/` inside the image. The latter is not writable when the container runs as a non-root UID (the standard configuration in this pipeline), and the upstream-staged library also contains a non-directory entry at `famdb`. The private libdir is populated with symlinks to the upstream library entries plus the Dfam partition files, and is passed to RepeatMasker via its `-libdir` argument.

## 9. Script migration

Several scripts from the original Snakemake workflow were adapted for execution inside Nextflow processes. Python scripts were migrated from Snakemake-specific inputs and outputs to command-line arguments using `sys.argv` or `argparse`. R scripts were migrated from the Snakemake API to `commandArgs` and explicit RMarkdown parameters where needed. The migration preserves analytical logic; only the I/O wiring was rewritten.

Adapted scripts include: `getGoodAlts.py`, `merge.py`, `genotype.py`, `vcfToBedForMosdepth.py`, `checkInsertionsMultiSample.py`, `getMeDeletions.py`, `svaf.py`, `analysisPreparatory.R`, `analysisGenotyping.R`, `enrichment.R`, `report.Rmd`, `comparison.Rmd`. A new `hallmarks.py` script was added for the retrotransposition hallmark detection module.

## 10. Main differences from the original pipeline

RetroInspector-NF introduces the following changes relative to the original RetroInspector workflow:

- Migration from Snakemake to Nextflow DSL2, with modular processes.
- Reference genome support extended from hg38 / GRCh38 only to hg38 / GRCh38 and T2T-CHM13v2.0.
- CHM13 RepeatMasker BED support and T2T-specific RepeatMasker library preparation.
- Adapted parsing of hg38 and CHM13 RepeatMasker BED formats where required.
- T2T-specific VCF header derived from CHM13v2.0 contig lengths.
- Dfam download anchored to release 3.9 to guarantee FamDB schema compatibility with the bundled RepeatMasker version, regardless of upstream alias changes.
- Optional retrotransposition hallmark detection, including polyA/polyT tails and L1 endonuclease motif detection.
- T2T-specific reporting components.

## 11. Troubleshooting

A few common issues seen during validation:

- **Docker image pull fails.** Check outbound connectivity to Docker Hub and confirm that `docker pull monicazh/retroinspector-nf:0.1.0` works outside Nextflow. On shared HPC nodes, Docker may be restricted; in that case use the `conda` profile or request a node with Docker enabled.
- **FTP download from 1000 Genomes is slow or interrupted.** The included scripts use `wget -c` so partial downloads are resumed automatically. If the FTP server is unreachable, mirror the tarballs locally and edit the script's `BASE_URL` variable.
- **`corruption of N bytes: Invalid chunk length` at Nextflow startup.** This is a glibc warning from a subprocess, usually triggered by mixing two conda/mamba installations on the host PATH. It does not affect the pipeline, but to silence it make sure the activated environment and the `mamba` on `PATH` come from the same installer.
- **Out-of-memory on REPEATMASKER or MINIMAP2.** Lower `cpus` and `memory` in the corresponding `withName` block of `nextflow.config`.

## Citation

This work is based on the original RetroInspector pipeline:

Cuenca-Guardiola J, de la Morena-Barrio B, Corral J, Fernández-Breis JT. *Advanced analysis of retrotransposon variation in the human genome with nanopore sequencing using RetroInspector.* Scientific Reports. 2025;15:14489. https://doi.org/10.1038/s41598-025-98847-7

---

Mónica Zapata Hernández — UMU
