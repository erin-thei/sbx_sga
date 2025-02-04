<img src="https://github.com/sunbeam-labs/sunbeam/blob/main/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_sga

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_sga/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_sga/actions/workflows/tests.yml)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_sga)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_sga/)
<!-- Badges end -->

## Introduction

sbx_sga (Single Genome Assembly) is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for isolate QC, assembly, and classification. This pipeline uses [Mash](https://github.com/marbl/mash) for quality control, [Shovill](https://github.com/tseemann/shovill) for bacterial isolate assembly, [CheckM2](https://github.com/chklovski/CheckM2) and [QUAST](https://github.com/ablab/quast) for assembly QC, [MLST](https://github.com/tseemann/mlst) for typing, [Bakta](https://github.com/oschwengers/bakta) for annotation, [abriTAMR](https://github.com/MDU-PHL/abritamr) for AMR profiling, and [Sylph](https://github.com/bluenote-1577/sylph) for taxonomic classification.

## Installation

Extension install is as simple as passing the extension's URL on GitHub to `sunbeam extend`:

    sunbeam extend https://github.com/sunbeam-labs/sbx_sga

Any user-modifiable parameters specified in `config.yml` are automatically added on `sunbeam init`. If you're installing an extension in a project where you already have a config file, run the following to add the options for your newly added extension to your config (the `-i` flag means in-place config file modification; remove the `-i` flag to see the new config in stdout):

    sunbeam config update -i /path/to/project/sunbeam_config.yml

Installation instructions for older versions of Sunbeam are included at the end of this README.

## Running

To run an extension, simply run Sunbeam as usual with your extension's target rule specified:

    sunbeam run --profile /path/to/project/ example_rule

### Options for config.yml

  - mash_ref: the reference file for running Mash (should be a file ending in `.msh`)
  - checkm_ref: the diamond database for running CheckM2 (should be a file ending in `.dmnd`)
  - bakta_ref: the bakta reference database (should be a directory similar to `.../bakta_db/db/`)
    
## Installing an extension (legacy instructions for sunbeam <3.0)

Installing an extension is as simple as cloning (or moving) your extension directory into the sunbeam/extensions/ folder, installing requirements through Conda, and adding the new options to your existing configuration file: 

    git clone https://github.com/sunbeam-labs/sbx_sga/ sunbeam/extensions/sbx_sga
    cat sunbeam/extensions/sbx_sga/config.yml >> sunbeam_config.yml

## Issues with pipeline

Please post any issues with this extension [here](https://github.com/sunbeam-labs/sbx_sga/issues).
