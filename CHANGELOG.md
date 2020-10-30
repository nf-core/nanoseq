# nf-core/nanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

* Transcript reconstruction and quantification ([`bambu`](https://github.com/GoekeLab/bambu) or [`StringTie2`](https://ccb.jhu.edu/software/stringtie/) and [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
* Differential expression analysis at the gene-level ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and transcript-level ([`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html))
* Ability to provide BAM input to the pipeline
* Change samplesheet format to be more flexible to BAM input files
* Add AWS full-sized test data
* Add parameter JSON schema for pipeline
* Add citations file
* Add pycoQC and featureCounts output to MultiQC report
* Update pipeline template to nf-core/tools `1.11`
* Collapsible sections for output files in `docs/output.md`
* Replace `set` with `tuple` and `file()` with `path()` in all processes
* Capitalise process names
* Parameters:
  * `--barcode_both_ends` Require barcode on both ends for Guppy basecaller
  * `--quantification_method` to specify the transcript quantification method to use
  * `--skip_quantification` to skip transcript quantification and differential analysis
  * `--skip_differential_analysis` to skip differential analysis with DESeq2 and DEXSeq
  * `--publish_dir_mode` to customise method of publishing results to output directory [nf-core/tools#585](https://github.com/nf-core/tools/issues/585)

### `Fixed`

* Added `--gpus all` to Docker `runOptions` when using GPU as mentioned [here](https://github.com/docker/compose/issues/6691#issuecomment-514429646)

### `Dependencies`

* Add seaborn `0.10.1`
* Add r-base `4.0.3`
* Add stringtie `2.1.4`
* Add subread `2.0.1`
* Add bioconductor-deseq2 `1.28.0`
* Add bioconductor-drimseq `1.16.0`
* Add bioconductor-dexseq `1.34.0`
* Add bioconductor-stager `1.10.0`
* Add bioconductor-bsgenome `1.56.0`
* Update python `3.7.3` -> `3.8.6`
* Update markdown `3.1.1` -> `3.3.3`
* Update pymdown-extensions `6.0` -> `8.0.1`
* Update pygments `2.5.2` -> `2.7.2`
* Update multiqc `1.8` -> `1.9`
* Update nanoplot `1.28.4` -> `1.32.1`
* Update samtools `1.9` -> `1.11`
* Update ucsc-bedgraphtobigwig `357` -> `377`
* Update ucsc-bedtobigbed `357` -> `377`
* Remove psutil `5.7.0`

## [1.0.0] - 2020-03-05

Initial release of nf-core/nanoseq, created with the [nf-core](http://nf-co.re/) template.
