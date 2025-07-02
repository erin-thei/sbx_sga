# Test Data

This directory contains read files and miniature reference databases used for testing `sbx_sga`.

## Sample reads

The `reads/` subdirectory holds paired-end FASTQ files generated from two bacterial genomes with a handful of host reads. The genomes are available from the [Sunbeam test data repository](https://github.com/sunbeam-labs/sunbeam/tree/main/tests/data/raw). The reads were produced using `art_illumina` with 150 bp read length.

## Miniature databases

The end-to-end tests can be run with very small versions of the databases required by the pipeline. Each database should be compact so that it can be stored under version control and executed on GitHub Actions.

## Setup

Create a small conda environment containing the tools needed to build the test
databases:

```bash
conda create -n sga_testdbs -c conda-forge -c bioconda \
    mash bakta checkm2 genomad sylph diamond prodigal
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
   ```bash
   prodigal -i genome1.fna -a genome1.faa
   prodigal -i genome2.fna -a genome2.faa
   cat genome1.faa genome2.faa > species.faa
   ```
2. Create a DIAMOND database:
   ```bash
   diamond makedb --in species.faa -d checkm_test.dmnd
   ```

### Bakta
Bakta has a number of different database components that each present their own challenges for miniaturization.

### Genomad
Genomad is another tricky one. Haven't found a way yet to create minimal versions of all the required components.

### Sylph
1. Construct a minimal Sylph database from the two genomes:
   ```bash
   sylph sketch genome1.fna genome2.fna
   ```
2. This produces a lightweight index suitable for `sylph classify`, keeping the directory under a few megabytes.

Place the resulting files in `.tests/data/` and configure the tests to point `mash_ref`, `checkm_ref`, `bakta_ref`, `genomad_ref`, and `sylph_ref` to these paths.
