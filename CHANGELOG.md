# nf-core/nanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unpublished Version / DEV]

### `Added`

### `Fixed`

* Added `--gpus all` to Docker `runOptions` when using GPU as mentioned [here](https://github.com/docker/compose/issues/6691#issuecomment-514429646)

### `Dependencies`

* Add r-base `4.0.1`
* Add r-devtools `2.3.0`
* Add stringtie `2.0`
* Add subread `2.0.1`
* Add bioconductor-deseq2 `1.28.0`
* Add bioconductor-drimseq `1.16.0`
* Add bioconductor-dexseq `1.34.0`
* Add bioconductor-stager `1.10.0`
* Add bioconductor-bsgenome `1.56.0`

## [1.0.0] - 2020-03-05

Initial release of nf-core/nanoseq, created with the [nf-core](http://nf-co.re/) template.
