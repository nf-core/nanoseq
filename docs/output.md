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

