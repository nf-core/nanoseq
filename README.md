# ![nfcore/nanoseq](docs/images/nf-core-nanoseq_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/nanoseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/nanoseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/nanoseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/nanoseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/nanoseq.svg)](https://hub.docker.com/r/nfcore/nanoseq)

## Introduction

**nfcore/nanoseq** is a bioinformatics analysis pipeline that can be used to perform basecalling, demultiplexing, mapping and QC of Nanopore DNA/RNA sequencing data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary

1. Basecalling and/or demultiplexing ([`Guppy`](https://nanoporetech.com/nanopore-sequencing-data-analysis) or [`qcat`](https://github.com/nanoporetech/qcat); *optional*)
2. Sequencing QC ([`pycoQC`](https://github.com/a-slide/pycoQC), [`NanoPlot`](https://github.com/wdecoster/NanoPlot))
3. Raw read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot), [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Alignment ([`GraphMap2`](https://github.com/lbcb-sci/graphmap2) or [`minimap2`](https://github.com/lh3/minimap2))
    * Both aligners are capable of performing unspliced and spliced alignment. Sensible defaults will be applied automatically based on a combination of the input data and user-specified parameters
    * Each sample can be mapped to its own reference genome if multiplexed in this way
    * Convert SAM to co-ordinate sorted BAM and obtain mapping metrics ([`SAMtools`](http://www.htslib.org/doc/samtools.html))
5. Create bigWig ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/)) and bigBed ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedToBigBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)) coverage tracks for visualisation
6. Present QC for alignment results ([`MultiQC`](https://multiqc.info/docs/))
7. Transcript reconstruction and quantification ([`bambu`](https://github.com/GoekeLab/bambu) or [`StringTie2`](https://ccb.jhu.edu/software/stringtie/))
   * When [`StringTie2`](https://ccb.jhu.edu/software/stringtie/) is chosen, each sample can be processed individually and combined. After which, [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/) will be used for both gene and transcript quantification.
   * [`bambu`](https://github.com/GoekeLab/bambu) performs both transcript reconstruction and quantification.
8. Differential expression analysis with ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) for condition comparison
   * At least 3 replicates for each condtion need to be satistified for this step. 

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/) or [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/nanoseq -profile test,<docker/singularity/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/nanoseq \
    --input samplesheet.csv \
    --protocol DNA \
    --input_path ./fast5/ \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --barcode_kit SQK-PBK004 \
    -profile <docker/singularity/institute>
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline. An example input samplesheet for performing both basecalling and demultiplexing can be found [here](assets/samplesheet.csv).

## Documentation

The nf-core/nanoseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

nf-core/nanoseq was originally written by [Chelsea Sawyer](https://github.com/csawye01) and [Harshil Patel](https://github.com/drpatelh) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) for use at [The Francis Crick Institute](https://www.crick.ac.uk/), London. Other primary contributors include [Laura Wratten](https://github.com/lwratten), [Ying Chen](https://github.com/cying111), [Yuk Kei Wan](https://github.com/yuukiiwa) and [Jonathan Goeke](https://github.com/jonathangoeke) from the [Genome Institute of Singapore](https://www.a-star.edu.sg/gis), [Johannes Alneberg](https://github.com/alneberg) and [Franziska Bonath](https://github.com/FranBonath) from [SciLifeLab](https://www.scilifelab.se/), Sweden.

Many thanks to others who have helped out along the way too, including (but not limited to): [@crickbabs](https://github.com/crickbabs), [@AnnaSyme](https://github.com/AnnaSyme).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/nanoseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use  nf-core/nanoseq for your analysis, please cite it using the following doi: [10.5281/zenodo.3697959](https://doi.org/10.5281/zenodo.3697959)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
