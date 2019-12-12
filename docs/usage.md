# nf-core/nanoseq: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--input`](#--input)
  * [`--protocol`](#--protocol)
* [Basecalling](#basecalling)
  * [`--run_dir`](#--run_dir)
  * [`--flowcell`](#--flowcell)
  * [`--kit`](#--kit)
  * [`--barcode_kit`](#--barcode_kit)
  * [`--guppy_config`](#--guppy_config)
  * [`--guppy_model`](#--guppy_model)
  * [`--guppy_gpu`](#--guppy_gpu)
  * [`--guppy_gpu_runners`](#--guppy_gpu_runners)
  * [`--guppy_cpu_threads`](#--guppy_cpu_threads)
  * [`--gpu_device`](#--gpu_device)
  * [`--gpu_cluster_options`](#--gpu_cluster_options)
  * [`--skip_basecalling`](#--skip_basecalling)
  * [`--skip_demultiplexing`](#--skip_demultiplexing)
* [Alignments](#alignments)
  * [`--stranded`](#--stranded)
  * [`--aligner`](#--aligner)
  * [`--save_align_intermeds`](#--save_align_intermeds)
  * [`--skip_alignment`](#--skip_alignment)
* [Coverage tracks](#coverage-tracks)
* [Skipping QC steps](#skipping-qc-steps)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)
<!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

A typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/nanoseq \
    --input samplesheet.csv \
    --protocol DNA \
    --run_dir ./fast5/ \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --barcode_kit SQK-PBK004 \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/nanoseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/nanoseq releases page](https://github.com/nf-core/nanoseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/nanoseq`](http://hub.docker.com/r/nfcore/nanoseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/nanoseq`](http://hub.docker.com/r/nfcore/nanoseq/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--input`

You will need to create a file with information about the samples in your experiment/run before executing the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row. As shown in the examples below, the accepted format of the file is slightly different if you would like to run the pipeline with or without demultiplexing.

#### With basecalling and demultiplexing

```bash
sample,fastq,barcode,genome
Sample1,,1,mm10
Sample2,,2,mm10
Sample3,,3,hg19
Sample4,,4,/path/to/local/reference/genome.fa
```

#### With basecalling but not demultiplexing

```bash
sample,fastq,barcode,genome
Sample1,,1,mm10
```

> You will have to specify the `--skip_demultiplexing` parameter if you wish to bypass the demultiplexing step.

#### Without both basecalling and demultiplexing

```bash
sample,fastq,barcode,genome
Sample1,SAM101A1.fastq.gz,,mm10
Sample2,SAM101A2.fastq.gz,,mm10
Sample3,SAM101A3.fastq.gz,,hg19
Sample4,SAM101A4.fastq.gz,,/path/to/local/reference/genome.fa

> You will have to specify the `--skip_basecalling` parameter if you wish to bypass the basecalling and demultiplexing steps.

```

| Column   | Description                                                                                                                |
|----------|----------------------------------------------------------------------------------------------------------------------------|
| `sample` | Sample name without spaces.                                                                                                |
| `fastq`  | Full path to FastQ file if previously demultiplexed. File has to be zipped and have the extension ".fastq.gz" or ".fq.gz". |
| `barcode`| Barcode identifier attributed to that sample when multiplexing samples in integer format.                                  |
| `genome` | Genome fasta for alignment. This can either be a local path, or the appropriate key for a genome available on [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) (see [iGenomes config file](../conf/igenomes.config)). If unspecified then the alignment step will be skipped for that sample. |

### `--protocol`

Specifies the type of data that was sequenced i.e. "DNA", "cDNA" or "directRNA".

## Basecalling

### `--run_dir`

Path to Nanopore run directory e.g. `fastq_pass/`

### `--flowcell`

Flowcell used to perform the sequencing e.g. `FLO-MIN106`. Not required if `--guppy_config` is specified.

### `--kit`

Kit used to perform the sequencing e.g. `SQK-LSK109`. Not required if `--guppy_config` is specified.

### `--barcode_kit`

Barcode kit used to perform the sequencing e.g. `SQK-PBK004`

### `--guppy_config`

Guppy config file used for basecalling passed with the `--config` parameter. Cannot be used in conjunction with `--flowcell` and `--kit`.
This can be a local file (i.e. `/your/dir/guppy_conf.cfg`) or a string specifying a configuration stored in the `/opt/ont/guppy/data` directory of Guppy.

### `--guppy_model`

Custom basecalling model file (`json`) to pass to Guppy for basecalling with the `--model` parameter. Custom basecalling models can be trained with software such as [Taiyaki](https://github.com/nanoporetech/taiyaki). This can also be a string specifying a model stored in the `/opt/ont/guppy/data` directory of Guppy.

### `--guppy_gpu`

Whether to demultiplex with Guppy in GPU mode.

### `--guppy_gpu_runners`

Number of '--gpu_runners_per_device' used for guppy when using `--guppy_gpu` (default: 6)

### `--guppy_cpu_threads`

Number of '--cpu_threads_per_caller' used for guppy when using `--guppy_gpu` (default: 1)

### `--gpu_device`

Basecalling device specified to Guppy in GPU mode using `--device` (default: 'auto')

### `--gpu_cluster_options`

Cluster options required to use GPU resources (e.g. '--part=gpu --gres=gpu:1')

### `--skip_basecalling`

Skip basecalling with Guppy

### `--skip_demultiplexing`

Skip demultiplexing with Guppy

## Alignment

### `--stranded`

Specifies if the data is strand-specific. Automatically activated when using --protocol directRNA (default: false)

When using `--protocol`/`--stranded` the following command-line arguments will be set for `minimap2` and `graphmap2`:

| `nanoseq` input              | `minimap2` presets  | `graphmap2` presets |
|------------------------------|---------------------|--------------------|
| `--protocol DNA`             | -ax map-ont         | no presets         |
| `--protocol cDNA`            | -ax splice          | -x rnaseq          |
| `--protocol directRNA`       | -ax splice -uf -k14 | -x rnaseq          |
| `--protocol cDNA --stranded` | -ax splice -uf      | -x rnaseq          |

### `--aligner`

Specifies the aligner to use (available are: `graphmap2` or `minimap2`)

### `--save_align_intermeds`

Save the `.sam` files from the alignment step - not done by default

### `--skip_alignment`

Skip alignment and subsequent process

## Coverage tracks

| Step                    | Description                            |
|-------------------------|----------------------------------------|
| `--skip_bigwig`         | Skip BigWig file generation            |
| `--skip_bigbed`         | Skip BigBed file generation            |

## Skipping QC steps

The pipeline contains a number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

| Step                    | Description                          |
|-------------------------|--------------------------------------|
| `--skip_qc`             | Skip all QC steps apart from MultiQC |
| `--skip_pycoqc`         | Skip pycoQC                          |
| `--skip_nanoplot`       | Skip NanoPlot                        |
| `--skip_fastqc`         | Skip FastQC                          |
| `--skip_multiqc`        | Skip MultiQC                         |

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
