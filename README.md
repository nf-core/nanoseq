# ![nfcore/nanodemux](docs/images/nfcore-nanodemux_logo.png)

[![Build Status](https://travis-ci.com/nf-core/nanodemux.svg?branch=master)](https://travis-ci.com/nf-core/nanodemux)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nanodemux.svg)](https://hub.docker.com/r/nfcore/nanodemux)

## Introduction

**nfcore/nanodemux** is a bioinformatics analysis pipeline that can be used to demultiplex, QC and map Nanopore data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary
This pipeline
1. Checking the sample sheet
    * Checks to see if the sample sheet is formatted correctly.
2. Basecalling and Barcoding 
    * [Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) if fastq files are not provided, this process performs demultiplexing by basecalling and barcoding.
    * Must have access to the Nanopore software Guppy and have the kit and barcode kit used on the sample sheet in the format of the example sample sheet provided.
    * Can be bypassed using the `--skipDemultiplexing` parameter and the inclusion of fastq file links in the sample sheet.
3. PycoQC and NanoPlot for quality checking the data
    * [PycoQC](https://github.com/a-slide/pycoQC) outputs diagnostic plots and data for quality control of sequencing data from [Oxford Nanopore Technology](https://nanoporetech.com/) in an HTML format.
    * [NanoPlot](https://github.com/wdecoster/NanoPlot) plotting tool for long read sequencing data and alignments.
    * Can be bypassed using the `--skipPycoQC`, `--skipNanoPlot` or `--skipQC` (for both PycoQC and NanoPlot) parameter.
4. Alignment (choice of either using Graphmap or Minimap2)
    * Choose the aligner with the parameter `--aligner`.
	* [GraphMap](https://github.com/isovic/graphmap) is a highly sensitive and accurate mapper for long, error-prone reads.
    * [Minimap2](https://github.com/lh3/minimap2) a versatile pairwise aligner for genomic and spliced nucleotide sequences.
    * Can be bypassed using the `--skipAlignment` parameter.
5. Sorting BAM files
	* [SAMtools](http://www.htslib.org/doc/samtools.html) Converting .bam files to coordinate sorted and indexed .bam.
6. MultiQC
    * [MultiQC](https://multiqc.info/docs/) Presenting the QC for the output from SAMtools: Percent mapped, Alignment metrics, Samtools flagstat, and Mapped reads per contig.  


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
