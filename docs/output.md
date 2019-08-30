# nf-core/nanodemux: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Guppy](#guppy) - demultiplexing of Nanopore data
* [PycoQC](#pycoqc) - read quality control
* [NanoPlot](#nanoplot) - read quality control
* [GraphMap](#graphmap) - mapping for long reads
* [MiniMap2](#minimap2) - mapping for long reads
* [SortBam](#sortbam) - coordinate sort BAM files using SAMtools
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline


## Demultiplexing
*Documentation*: 
[Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)

**Output directory: `results/guppy`**

## Quality Control 
*Documentation*:
[PycoQC](https://github.com/a-slide/pycoQC)
[NanoPlot](https://github.com/wdecoster/NanoPlot)

*Description*:


*Output directories*: 
* `results/pycoQC`
pycoQC_output.html
* `results/nanoplot/fastq`

## Alignment 
*Documentation*:
[GraphMap](https://github.com/isovic/graphmap)
[MiniMap2](https://github.com/lh3/minimap2)
[SortBam](http://www.htslib.org/doc/samtools.html)

*Output directories*:
 * `results/graphmap`
 * `results/minimap2`

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
