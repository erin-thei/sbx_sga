# Test Data

This directory contains read files and miniature reference databases used for testing `sbx_sga`.

## Sample reads

The `reads/` subdirectory holds paired-end FASTQ files generated from two bacterial genomes with a handful of host reads. The genomes are available from the [Sunbeam test data repository](https://github.com/sunbeam-labs/sunbeam/tree/main/tests/data/raw). The reads were produced using `art_illumina` with 150&nbsp;bp read length.

## Miniature databases

The end-to-end tests can be run with very small versions of the databases required by the pipeline. Each database should be compact so that it can be stored under version control and executed on GitHub Actions.

## Setup

Create a small conda environment containing the tools needed to build the test
databases:

```bash
conda create -n sga_testdbs -c conda-forge -c bioconda \
    mash bakta checkm2 genomad diamond prodigal
conda activate sga_testdbs
```

### Mash
1. Download the two bacterial genomes used to create the reads.
2. Build a small sketch:
   ```bash
   mash sketch -o mash_test.msh genome1.fna genome2.fna
   ```

### CheckM2
1. Extract protein sequences from the genomes (e.g. with `prodigal`).
2. Create a DIAMOND database:
   ```bash
   diamond makedb --in proteins.faa -d checkm_test.dmnd
   ```

### Bakta
1. Create a new Bakta database populated only with the proteins from the test genomes. This keeps the data set tiny:
   ```bash
   mkdir bakta_db
   bakta database init bakta_db
   bakta database add-genomes genome1.fna genome2.fna --output bakta_db
   ```
2. Only the core taxonomy files and the two genomes are included so the directory remains under a few megabytes.

### Genomad
1. Build a miniature Genomad database using the same bacterial references:
   ```bash
   mkdir genomad_db
   genomad build-database genome1.fna genome2.fna genomad_db/
   ```
2. The command generates only the indices required for `genomad end-to-end`, producing a very small database directory.

Place the resulting files in a directory such as `.tests/data/databases/` and configure the tests to point `mash_ref`, `checkm_ref`, `bakta_ref`, and `genomad_ref` to these paths.
