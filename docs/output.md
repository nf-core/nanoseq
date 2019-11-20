# nf-core/nanoseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Guppy](#guppy) - demultiplexing of Nanopore data
* [PycoQC](#pycoqc) - read quality control
* [NanoPlot](#nanoplot) - read quality control
* [GraphMap](#graphmap) - mapping for long reads
* [MiniMap2](#minimap2) - mapping for long reads
* [SortBam](#sortbam) - coordinate sort BAM files using SAMtools
* [bedtools](#bedtools) - create bigWig and bigBed files 
* [MultiQC](#multiqc) - aggregate report, describing results of the alignment

## Demultiplexing

*Documentation*:  
[Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)

*Description*:  
Guppy will demultiplex and barcode the data given from an ONT device. The flowcell, kit and barcode kit must be given in the command line if demultiplexing needed. This step can by bypassed using the `--skip_demultiplexing` parameter when initiating the pipeline. The output folders will be separated into the barcodes from the kit used and unclassified. The output in each barcode folder is then merged into one fastq file for easier downstream processing.

*Output directories*:

* `guppy/basecalling/barcode*/`  
  FastQ files output for each barcode
* `guppy/basecalling/unclassified/`  
  FastQ files output that are unclassified
* `guppy/fastq`
  Merged output of fastq files into one fastq for each barcode

## Sequencing Quality Control

*Documentation*:  
[PycoQC](https://github.com/a-slide/pycoQC), [NanoPlot](https://github.com/wdecoster/NanoPlot)

*Description*:  
PycoQC and NanoPlot give general quality metrics about the sequencing run. It provides information about the distribution of read length, read length over time, number of reads per barcode and other general stats.
![PycoQC - Number of Reads per Barcode plot](images/NumberofReadsperBarcode.png)

*Output directories*:

* `pycoQC/`  
  An .html file output is produced that includes a run summary and graphical representation of distribution of read length, distribution of read quality scores, mean read quality per sequence length, output per channel over experiment time, output over experiment time, read quality over experiment time, readlength over experiment time, and percentage of reads per barcode.
* `nanoplot/summary/`  
  An output of a .png file of a statistical run summary and an html summary file of overall run.


## FastQ Quality Control
*Documentation*:  
[NanoPlot](https://github.com/wdecoster/NanoPlot)

*Description*:  
NanoPlot give general quality metrics about the fastq output per barcode from Guppy. It provides information about the quality score distribution across your reads, read lengths and other general stats.
![Nanoplot - Read quality vs read length](images/NanoPlot_output.png)

*Output directories*:

* `nanoplot/fastq/`  
  An output of QC metric plots in individual .png files and in one html file summarizing the output. 

## Alignment

*Documentation*:  
[GraphMap](https://github.com/isovic/graphmap), [MiniMap2](https://github.com/lh3/minimap2), [SortBam](http://www.htslib.org/doc/samtools.html)

*Description*:  
The FastQ reads are mapped to the given reference assembly provided using either GraphMap or Minimap2 and then sorted and indexed using SAMtools or these processes can be bypassed using the `--skip_alignment` parameter.

The files resulting from the alignment with graphmap or minimap2 of individual libraries are not saved by default so this directory will not be present in your results. You can override this behaviour with the use of the `--save_align_intermeds` flag in which case it will contain the coordinate sorted alignment files in [`*.bam`](https://samtools.github.io/hts-specs/SAMv1.pdf) format.

![ALIGNER - Alignment per barcode](images/mqc_samtools_alignment_plot_1.png)

*Output directories*:

* `graphmap/`  
  If the `--aligner graphmap` parameter is used, the sorted and indexed bam files will be output here.
* `minimap2/`  
  If the `--aligner minimap2` parameter is used, the sorted and indexed bam files will be output here.
* `<ALIGNER>/samtools_stats/`  
  `*.flagstat`, `*.idxstats` and `*.stats` files generated from the alignment files using SAMtools.

## bigWig and bigBed 

*Documentation*:
[`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/), [`bedToBigBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)

*Description*:
Creation of bigWig and bigBed coverage tracks for visualisation. This can be bypassed by setting the parameters `--skip_bigwig` and/or `--skip_bigbed`.

*Output directories*:
* `<ALIGNER>/bigwig`  
  The bigWig files will be output here.
* `<ALIGNER>/bigbed`  
  The bigbed files will be output here.

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

*Output directories*:
  
* `multiqc/Project_multiqc_report.html`  
  MultiQC report - a standalone HTML file that can be viewed in your web browser
* `multiqc/multiqc_data/`  
  Directory containing parsed statistics from the different tools used in the pipeline
* `multiqc/multiqc_plots/`
  Directory containing the image files of the graphs included in MultiQC

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
