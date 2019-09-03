# ![nfcore/nanodemux](docs/images/nfcore-nanodemux_logo.png)

[![Build Status](https://travis-ci.com/nf-core/nanodemux.svg?branch=master)](https://travis-ci.com/nf-core/nanodemux)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nanodemux.svg)](https://hub.docker.com/r/nfcore/nanodemux)

## Introduction

**nfcore/nanodemux** is a bioinformatics analysis pipeline that can be used to demultiplex, QC and map Nanopore data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->
## Pipeline Summary

1. Basecalling, barcoding and adapter trimming ([`Guppy`](https://nanoporetech.com/nanopore-sequencing-data-analysis); *optional*)
    * Sequencing QC ([`pycoQC`](https://github.com/a-slide/pycoQC), [`NanoPlot`](https://github.com/wdecoster/NanoPlot))
2. FastQ QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
3. Alignment ([`GraphMap`](https://github.com/isovic/graphmap) or [`minimap2`](https://github.com/lh3/minimap2))
    * Convert `.sam` to co-ordinate sorted `.bam` and obtain mapping metrics ([`SAMtools`](http://www.htslib.org/doc/samtools.html))
4. Present QC for alignment results ([`MultiQC`](https://multiqc.info/docs/))

## Documentation

The nf-core/nanodemux pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits
nf-core/nanodemux was originally written by [Chelsea Sawyer](https://github.com/csawye01) and [Harshil Patel](https://github.com/drpatelh) from The Bioinformatics & Biostatistics Group for use at The Francis Crick Institute, London.

Many thanks to others who have helped out along the way too, including (but not limited to):  [`@crickbabs`](https://github.com/crickbabs)

## Citation

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  {{ cookiecutter.name }} for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
