# nf-core/nanodemux

**A pipeline to demultiplex ONT Nanopore data, quality control and simple alignment**.

[![Build Status](https://travis-ci.com/nf-core/nanodemux.svg?branch=master)](https://travis-ci.com/nf-core/nanodemux)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nanodemux.svg)](https://hub.docker.com/r/nfcore/nanodemux)

## Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.


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

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->
## Pipeline Summary
This pipeline 
1. Checking the sample sheet
    * Searches for [Data] tag and finds the column Fastq_path. If this column is empty it will send the run directory to the demultiplexing process and if not it bypasses to the QC step
2. Basecalling and Barcoding (CONDITIONAL)
    * If fastq files are not provided, this process performs demultiplexing by basecalling and barcoding.
    * Must have access to the Nanopore software Guppy and have the kit and barcode kit used on the sample sheet in the format of the example sample sheet provided
3. pycoQC for quality checking the data
    * [MinIONQC](https://github.com/roblanf/minion_qc) outputs diagnostic plots and data for quality control of sequencing data from Oxford Nanopore's MinION in an HTML format.
4. GraphMap
	* [GraphMap](https://github.com/isovic/graphmap) is a highly sensitive and accurate mapper for long, error-prone reads
5. Sorting BAM files 
	* Converting .bam files to coordinate sorted .bam
6. MultiQC
    * [MultiQC](https://multiqc.info/docs/)  

## Credits
nf-core/nanodemux was originally written by [Chelsea Sawyer](https://github.com/csawye01) and [Harshil Patel](https://github.com/drpatelh) from The Bioinformatics & Biostatistics Group for use at The Francis Crick Institute, London.

Many thanks to others who have helped out along the way too, including (but not limited to):  [`@crickbabs`](https://github.com/crickbabs)

## Citation

R Lanfear, M Schalamun, D Kainer, W Wang, B Schwessinger (2018). **MinIONQC: fast and simple quality control for MinION sequencing data**, Bioinformatics, [bty654](https://doi.org/10.1093/bioinformatics/bty654)

Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).