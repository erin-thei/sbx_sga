# Test Data

This directory contains read files and miniature reference databases used for testing `sbx_sga`.

## Sample reads

The `reads/` subdirectory holds paired-end FASTQ files generated from two bacterial genomes with a handful of host reads. The genomes are available from the [Sunbeam test data repository](https://github.com/sunbeam-labs/sunbeam/tree/main/tests/data/raw). The reads were produced using `art_illumina` with 150&nbsp;bp read length.

## Miniature databases

The end-to-end tests can be run with very small versions of the databases required by the pipeline. Each database should be compact so that it can be stored under version control and executed on GitHub Actions.

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
1. Download the Bakta database into a new directory:
   ```bash
   bakta database --download bakta_db/
   ```
2. Remove large optional files to keep the directory small (taxdump archives, large FASTA files, etc.).

### Genomad
1. Obtain the Genomad database:
   ```bash
   genomad download-database genomad_db/
   ```
2. Delete nonessential files to reduce its size, leaving `version.txt` and the necessary indices.

Place the resulting files in a directory such as `.tests/data/databases/` and configure the tests to point `mash_ref`, `checkm_ref`, `bakta_ref`, and `genomad_ref` to these paths.
