# FungiFunk

**Fungal functional profiling from paired-end sequencing reads.**

FungiFunk is a modular command-line toolkit for functional profiling of fungal paired-end sequencing reads.

It combines protein-level read mapping with DIAMOND, EggNOG-derived functional annotation, multi-sample processing, cross-sample summaries, and optional construction of custom reference databases.

Its modules can be used independently or connected: a custom reference can be built with `FungiFunk_build`, multiple samples can be processed with `FungiFunk_batch`, which runs the core `FungiFunk` workflow for each sample, and the resulting outputs can then be combined with `FungiFunk_summary`.

**Current version:** v0.1.9\
**Status:** Active development

## How FungiFunk fits together

FungiFunk is organised as a set of connected modules covering different stages of a functional profiling workflow.

A complete multi-sample workflow can use the modules as follows:

```text
                         Reference database
                    ┌────────────┴────────────┐
                    │                         │
       Official FungiFunk database*     Custom protein FASTA
                    |                          │
                    |                          ▼
                    |               FungiFunk_build.sh
                    |                          │
                    |               Paired DIAMOND + EggNOG
                    |                  reference databases
                    │                         │
                    └────────────┬────────────┘
                                 ▼
                       Paired-end FASTQ reads
                                 │
                      Multiple samples + TSV
                                 │
                                 ▼
                      FungiFunk_batch.sh
                                 │
                      runs FungiFunk.sh
                         for each sample
                                 │
                                 ▼
                       Per-sample outputs
                                 │
                                 ▼
                     FungiFunk_summary.sh
                                 │
                                 ▼
             Cross-sample count matrices
                  + summary tables
             + interactive HTML report
```

*An official FungiFunk reference database set is planned for distribution. Once available, users will also be able to skip the custom database-building step and use the provided reference databases directly.

The modules can also be used independently. For example, `FungiFunk.sh` can analyse a single sample directly, and `FungiFunk_summary.sh` can combine outputs from independently executed FungiFunk runs without requiring `FungiFunk_batch.sh`.

### Main modules

| Module | Purpose |
| --- | --- |
| `FungiFunk_build.sh` | Build a compatible custom DIAMOND + EggNOG reference database pair from a protein FASTA file. |
| `FungiFunk.sh` | Perform functional profiling for one paired-end sample. |
| `FungiFunk_batch.sh` | Run `FungiFunk.sh` sequentially across multiple samples from a TSV samplesheet. |
| `FungiFunk_summary.sh` | Combine outputs from multiple FungiFunk runs into cross-sample matrices, summary tables, and an interactive HTML report. |

The Bash scripts are the recommended entry points. Python scripts in the repository are used internally by these modules.

### What happens during sample profiling?

For each sample, `FungiFunk.sh`:

1. Takes paired-end sequencing reads in FASTQ or compressed FASTQ format.
2. Maps R1 and R2 reads independently against a fungal protein reference database using DIAMOND `blastx`.
3. Links detected reference proteins to their EggNOG-derived functional annotation.
4. Produces a per-read annotation table.
5. Generates read-count tables for:
   - Reference proteins
   - KEGG orthologs (KOs)
   - eggNOG orthologous groups (OGs)
   - CAZy families
   - PFAM domains

The sample-level profiling step can be represented as:

```text
Paired-end reads
   │
   ├── R1 ──► DIAMOND blastx ──┐
   │                           │
   └── R2 ──► DIAMOND blastx ──┤
                               ▼
                    Protein reference hits
                               │
                               ▼
                    EggNOG annotation transfer
                               │
                               ▼
                     Functional annotation
                               │
        ┌──────────────┬────────┼────────┬──────────────┐
        ▼              ▼        ▼        ▼              ▼
    Proteins       KEGG KOs  eggNOG OGs  CAZy         PFAM
        │              │        │        │              │
        └──────────────┴────────┴────────┴──────────────┘
                               │
                               ▼
                       Read-count tables
```

---

# Installation

FungiFunk runs from the command line using Bash, Python, DIAMOND, EggNOG-mapper, and supporting Python packages.

All FungiFunk modules use the same Conda environment.

A environment configuration file, `FungiFunk.yml`, is provided in the root directory of the repository.

## 1. Get FungiFunk

Clone or download the repository and move into its main directory:

```bash
git clone https://github.com/Rounge-lab/FungiFunk
cd FungiFunk
```

## 2. Create the FungiFunk Conda environment

Make sure Conda is installed and available in your terminal.

Create the environment using the provided YAML file:

```bash
conda env create -f FungiFunk.yml
```

Activate it:

```bash
conda activate FungiFunk
```

## 3. Reference databases

FungiFunk requires two paired reference files:

- A DIAMOND protein database (`.dmnd`)
- An EggNOG annotation table (`.emapper.annotations`)

Both files must represent the same underlying set of reference proteins so that DIAMOND hits can be matched correctly to their functional annotation.

Their paths can be provided explicitly with:

```text
-d, --diamond_db
-e, --eggnog_db
```

FungiFunk is designed to use an official default paired database set when these options are not provided.

> **Database availability:** The official FungiFunk default reference databases are planned for distribution but are not yet available as part of the current development release. Until they are provided, users should supply compatible DIAMOND and EggNOG database files explicitly, or place compatible files in the default FungiFunk database locations.

A compatible custom database pair can also be generated directly with `FungiFunk_build.sh`.

## 4. Check the installation

Check the version of each module:

```bash
bash FungiFunk.sh --version
bash FungiFunk_batch.sh --version
bash FungiFunk_summary.sh --version
bash FungiFunk_build.sh --version
```

Detailed command-line help is available with:

```bash
bash FungiFunk.sh --help
bash FungiFunk_batch.sh --help
bash FungiFunk_summary.sh --help
bash FungiFunk_build.sh --help
```

---

# Quick start

The FungiFunk modules can be used independently, but they can also be connected.

The example below shows a typical multi-sample analysis using a custom reference database.

## 1. Build a custom reference database

Starting from a protein database FASTA file:

```bash
bash FungiFunk_build.sh \
    -i proteins.faa \
    -o custom_database \
    -p my_fungal_db \
    -t 8
```

This creates a paired database set:

```text
custom_database/
├── diamond/
│   └── my_fungal_db.dmnd
└── eggnog/
    └── my_fungal_db.emapper.annotations
```

Both files originate from the same protein FASTA and can therefore be used together by FungiFunk.

If an official or previously prepared compatible FungiFunk database is already available, this step can be skipped.

## 2. Prepare a samplesheet

For multiple samples, create a tab-separated file with the following columns:

```text
prefix	reads1	reads2
sampleA	/path/sampleA_R1.fastq.gz	/path/sampleA_R2.fastq.gz
sampleB	/path/sampleB_R1.fastq.gz	/path/sampleB_R2.fastq.gz
```

## 3. Run FungiFunk across the samples

Pass the databases generated by `FungiFunk_build.sh` directly to `FungiFunk_batch.sh`:

```bash
bash FungiFunk_batch.sh \
    -s samplesheet.tsv \
    -o batch_results \
    -d custom_database/diamond/my_fungal_db.dmnd \
    -e custom_database/eggnog/my_fungal_db.emapper.annotations \
    -t 8
```

`FungiFunk_batch.sh` runs `FungiFunk.sh` once for each sample and stores the corresponding FungiFunk results in separate sample directories.

## 4. Summarise the dataset

The complete batch output can then be passed directly to `FungiFunk_summary.sh`:

```bash
bash FungiFunk_summary.sh \
    --input_dir batch_results \
    --outdir summary_results
```

The summary module searches the input directory recursively and generates cross-sample count matrices, summary tables, feature metadata, and an interactive HTML report.

The complete workflow is therefore:

```text
protein FASTA
     │
     ▼
FungiFunk_build
     │
     ▼
paired reference databases
     │
     ▼
FungiFunk_batch
     │
     ├──► FungiFunk → sampleA results
     ├──► FungiFunk → sampleB results
     └──► FungiFunk → sampleC results
                         │
                         ▼
                 FungiFunk_summary
                         │
                         ▼
                cross-sample results
```

## Running a single sample

`FungiFunk.sh` can also be used directly without the batch module:

```bash
bash FungiFunk.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -p sampleA \
    -o results/sampleA \
    -d custom_database/diamond/my_fungal_db.dmnd \
    -e custom_database/eggnog/my_fungal_db.emapper.annotations \
    -t 8
```

Here:

- `-1` provides the R1 reads.
- `-2` provides the R2 reads.
- `-p` defines the sample name used in output files.
- `-o` defines the sample output directory.
- `-d` and `-e` select the paired reference databases.
- `-t` controls the number of threads.

FASTQ and compressed FASTQ (`.fastq.gz`) input files are supported.

Outputs from multiple independently executed `FungiFunk.sh` runs can also be combined later with `FungiFunk_summary.sh`.

---

# Main FungiFunk outputs

A standard run produces a run-information log, input read counts, DIAMOND results, annotation tables, and functional read-count tables.

```text
sample_output/
│
├── FungiFunk_<prefix>_run_info.log
├── <prefix>_input_read_counts.tsv
│
├── diamond/
│   ├── <prefix>_R1.diamond.out
│   ├── <prefix>_R1.diamond.log
│   ├── <prefix>_R2.diamond.out
│   └── <prefix>_R2.diamond.log
│
└── annotation/
    ├── <prefix>_annotation.tsv
    ├── <prefix>_unannotated_hits.tsv
    ├── <prefix>_protein_read_counts.tsv
    ├── <prefix>_ko_read_counts.tsv
    ├── <prefix>_eggnog_ogs_read_counts.tsv
    ├── <prefix>_cazy_read_counts.tsv
    └── <prefix>_pfam_read_counts.tsv
```

## Annotation table

The standard annotation table contains the main DIAMOND alignment information together with functional annotation including:

- Reference protein
- DIAMOND E-value and bit score
- Read mate
- EggNOG seed ortholog
- eggNOG orthologous groups
- Functional description
- Preferred gene/protein name
- Enzyme Commission number
- KEGG ortholog
- KEGG pathway
- CAZy annotation
- PFAM domains

A larger extended annotation table can optionally be written using:

```text
--write_full_annotation
```

## Read-count tables

FungiFunk generates separate count tables for several functional levels.

These counts represent **read-level evidence**. R1 and R2 reads are processed independently and counts are not collapsed into paired-read fragments.

For eggNOG orthologous groups, FungiFunk counts the most specific valid OG reported for each annotated hit while retaining the complete EggNOG OG hierarchy in the annotation table.

Individual count tables can be disabled with the corresponding `--skip_*_read_counts` options, or all count tables can be disabled with:

```text
--skip_count_tables
```

Large annotation outputs can optionally be compressed with:

```text
--compress_outputs
```

Read-count summary tables remain uncompressed.

---

# Alignment filtering

FungiFunk can pass several filtering options directly to DIAMOND `blastx`:

| Option | Description |
| --- | --- |
| `--evalue` | Maximum E-value for reported alignments. |
| `--min-score` | Minimum bit score for reported alignments. |
| `--id` | Minimum sequence identity percentage. |
| `--query-cover` | Minimum query coverage percentage. |
| `--top` | Keep hits whose scores fall within the specified percentage of the best hit for a query. |

In DIAMOND, `--min-score` takes precedence over `--evalue`.

These filters determine which alignments reach the annotation stage.

For the complete list of options and current defaults, use:

```bash
bash FungiFunk.sh --help
```

---

# Running multiple samples with FungiFunk_batch

`FungiFunk_batch.sh` runs `FungiFunk.sh` sequentially across multiple paired-end samples.

`FungiFunk_batch.sh` uses the same `FungiFunk` Conda environment created from `FungiFunk.yml`.

## Samplesheet

The input samplesheet must be a tab-separated file containing:

```text
prefix	reads1	reads2
sampleA	/path/sampleA_R1.fastq.gz	/path/sampleA_R2.fastq.gz
sampleB	/path/sampleB_R1.fastq.gz	/path/sampleB_R2.fastq.gz
```

## Run the batch

```bash
bash FungiFunk_batch.sh \
    -s samplesheet.tsv \
    -o batch_results \
    -t 8
```

Custom databases can also be supplied:

```bash
bash FungiFunk_batch.sh \
    -s samplesheet.tsv \
    -o batch_results \
    -d databases/fungi.dmnd \
    -e databases/fungi.emapper.annotations \
    -t 8
```

Most FungiFunk output and DIAMOND filtering options can also be passed through the batch module.

By default, the batch continues if an individual sample fails.

To stop at the first failed sample, use:

```text
--stop_on_error
```

## Batch outputs

Each sample receives its own FungiFunk output directory.

The batch module additionally writes:

```text
FungiFunk_batch_run_info.log
FungiFunk_batch_summary.tsv
```

`FungiFunk_batch_summary.tsv` records the status and exit code of each sample together with the corresponding log-file paths.

For all batch options:

```bash
bash FungiFunk_batch.sh --help
```

---

# Summarising multiple samples

`FungiFunk_summary.sh` combines outputs from multiple FungiFunk runs.

It accepts results produced either by `FungiFunk_batch` or by independently executed FungiFunk runs.

`FungiFunk_summary.sh` uses the same `FungiFunk` Conda environment created from `FungiFunk.yml`.

The input directory is searched recursively for compatible FungiFunk outputs.

## Basic usage

```bash
bash FungiFunk_summary.sh \
    --input_dir batch_results \
    --outdir summary_results
```

By default, FungiFunk_summary creates:

- Cross-sample raw read-count matrices
- Sample-level summary tables
- Annotation summary tables
- Detected-feature summaries
- Feature metadata tables
- An interactive HTML report

## Interactive HTML report

The HTML report is generated by default and provides an overview of the dataset together with functional feature tables and plots.

To disable it:

```text
--skip_html_report
```

By default, the HTML report displays the top 50 functional features to remain responsive.

This can be changed with:

```text
--html_top_n <int>
```

or disabled entirely with:

```text
--html_no_top
```

Displaying all features may create a large HTML file and substantially increase memory usage.

## CPM-normalised matrices

Additional CPM-normalised matrices can be generated using:

```text
--normalise_cpm
```

CPM is calculated separately for each functional count table using the total count sum for that sample within that table.

It does **not** use the original number of FASTQ input reads.

For annotation fields that can contain multiple assignments, CPM therefore represents counts per million annotation assignments.

## PFAM and eggNOG context

Optional dataset-derived context can be added for PFAM and eggNOG OG features with:

```text
--calculate_pfam_eggnog_context
```

This scans the annotation tables and identifies commonly associated protein names and descriptions to make these feature tables easier to interpret.

This information is derived from the current FungiFunk dataset and should not be interpreted as official PFAM or EggNOG database metadata.

Because this option requires additional scanning of annotation tables, it may increase runtime.

## Summary outputs

The main output structure is:

```text
summary_results/
│
├── FungiFunk_summary_report.html
│
├── matrices/
│   ├── protein_read_count_matrix.tsv
│   ├── ko_read_count_matrix.tsv
│   ├── eggnog_ogs_read_count_matrix.tsv
│   ├── cazy_read_count_matrix.tsv
│   └── pfam_read_count_matrix.tsv
│
└── summaries/
    ├── detected_output_files.tsv
    ├── sample_summary.tsv
    ├── annotation_summary.tsv
    ├── detected_features_summary.tsv
    ├── protein_metadata.tsv
    │
    └── feature_metadata/
        ├── protein_feature_metadata.tsv
        ├── ko_metadata.tsv
        ├── cazy_metadata.tsv
        ├── pfam_metadata.tsv
        └── eggnog_ogs_metadata.tsv
```

When `--normalise_cpm` is enabled, additional CPM matrices are written to `matrices/`.

By default, missing or inconsistent files produce warnings where possible and the summary continues using the available data.

Use:

```text
--strict
```

to stop instead when expected files are missing or inconsistent.

For all summary options:

```bash
bash FungiFunk_summary.sh --help
```

---

# Building a custom reference database

`FungiFunk_build.sh` creates a paired DIAMOND database and EggNOG annotation table from the same protein FASTA file.

Keeping both outputs derived from the same FASTA ensures that protein identifiers found by DIAMOND can be matched correctly to the EggNOG annotation table.

## Build a custom database

```bash
bash FungiFunk_build.sh \
    -i proteins.faa \
    -o custom_database \
    -p my_fungal_db \
    -t 8
```

The main outputs are:

```text
custom_database/
│
├── diamond/
│   └── my_fungal_db.dmnd
│
├── eggnog/
│   └── my_fungal_db.emapper.annotations
│
├── FungiFunk_build_run_info.log
│
└── logs/
```

Unless an existing EggNOG-mapper data directory is supplied with:

```text
--eggnog_data_dir
```

the required EggNOG data are downloaded to:

```text
custom_database/eggnog_data/
```

The resulting database pair can then be used with FungiFunk:

```bash
bash FungiFunk.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -p sampleA \
    -o results/sampleA \
    -d custom_database/diamond/my_fungal_db.dmnd \
    -e custom_database/eggnog/my_fungal_db.emapper.annotations
```

Existing outputs can be replaced during database building with:

```text
--force
```

For all database-building options:

```bash
bash FungiFunk_build.sh --help
```

---

# Choosing which modules to use

FungiFunk does not require every module to be used in every analysis.

| Goal | Recommended module(s) |
| --- | --- |
| Analyse one paired-end sample | `FungiFunk.sh` |
| Analyse multiple paired-end samples using the same settings | `FungiFunk_batch.sh` |
| Combine outputs from multiple FungiFunk runs | `FungiFunk_summary.sh` |
| Build a compatible custom reference database | `FungiFunk_build.sh` |
| Run a complete custom-reference multi-sample workflow | `FungiFunk_build.sh` → `FungiFunk_batch.sh` → `FungiFunk_summary.sh` |

`FungiFunk_batch.sh` is a convenience layer around the core `FungiFunk.sh`: it runs the same sample-level analysis sequentially for each sample listed in the samplesheet.

`FungiFunk_summary.sh` does not require the samples to have been processed through the batch module. Outputs from independently executed `FungiFunk.sh` runs can also be placed under a common input directory and summarised together.

---

# Reproducibility

FungiFunk writes run-information logs containing the software version, command, input files, database paths, and relevant parameters.

These logs are intended to make analysis settings easier to track and reproduce.

For reproducible analyses, record:

- The FungiFunk version
- The reference database version
- The command and filtering options used
- The software environment
- Any custom reference-database construction settings

---

# Important notes

- FungiFunk reports **read-level counts**, not reconstructed fragment-level counts.
- R1 and R2 reads are processed independently.
- DIAMOND filtering is performed before functional annotation.
- DIAMOND and EggNOG reference files used together must represent the same protein reference set.
- Functional annotations are transferred from the EggNOG annotation associated with each detected reference protein.
- CPM values generated by `FungiFunk_summary` are normalised using the count totals of each functional table, not the original FASTQ read count.
- FungiFunk is under active development. Output formats and functionality may evolve between development versions.

---

# Getting help

Each module provides detailed command-line documentation:

```bash
bash FungiFunk.sh --help
bash FungiFunk_batch.sh --help
bash FungiFunk_summary.sh --help
bash FungiFunk_build.sh --help
```

Version information can be displayed with the corresponding `--version` option.

---

# Development and maintenance

FungiFunk was originally developed by **Blanca Candelas Galavis** as part of her PhD research in the [Rounge group](https://www.mn.uio.no/bils/english/groups/rounge-group/index.html) at the University of Oslo (UiO), Norway.

The tool is currently maintained by **Blanca Candelas Galavis**.

Additional contributors and future maintainers will be acknowledged here as the project evolves.

## License

FungiFunk is distributed under the **GNU General Public License, version 3 (GPLv3)**.

See the [LICENSE](LICENSE) file for the full license text.

## Citation

Formal citation information will be added as the project develops.

In the meantime, please include the software name and version when reporting its use in scientific work.