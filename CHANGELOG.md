# nf-core/nanoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.0] - 2023-03-10

### Major enhancements

- Removed the `guppy` basecaller as distributing it via a docker image is a breach to EULA
- Bump minimum Nextflow version from 21.10.3 -> 21.10.3
- Update pipeline template to nf-core/tools `2.7.2`
- Update `bambu` version from `1.0.2` to `2.0.0`

### Parameters

- Removed `--flowcell` as `nanoseq` no longer supports basecalling
- Removed `--kit` as `nanoseq` no longer supports basecalling
- Removed `--guppy_config` as `nanoseq` no longer supports basecalling
- Removed `--guppy_model` as `nanoseq` no longer supports basecalling
- Removed `--guppy_gpu` as `nanoseq` no longer supports basecalling
- Removed `--guppy_gpu_runners` as `nanoseq` no longer supports basecalling
- Removed `--guppy_cpu_threads` as `nanoseq` no longer supports basecalling
- Removed `--output_demultiplex_fast5` as `nanoseq` no longer supports basecalling
- Removed `--skip_basecalling` as `nanoseq` no longer supports basecalling
- Removed `--skip_pycoqc` as `nanoseq` no longer supports basecalling

### Software dependencies

| Dependency           | Old version | New version |
| -------------------- | ----------- | ----------- |
| `bioconductor-bambu` | 2.0.0       | 3.0.8       |

## [3.0.0] - 2022-06-21

### Major enhancements

- Add DNA variant calling functionality
- Add RNA modification and fusion detection functionalities
- Add `demux_fast5` module to output demultiplexed fast5 files when `--output_demultiplex_fast5` is set
- Add `--trim_barcodes` in Guppy basecaller to trim the barcodes from output fastq
- Port pipeline to the updated Nextflow DSL2 syntax adopted on nf-core/modules
- Removed `--publish_dir_mode` as it is no longer required for the new syntax
- Bump minimum Nextflow version from 21.04.0 -> 21.10.3
- Update pipeline template to nf-core/tools `2.2`
- Update `bambu` version from `1.0.2` to `2.0.0`
- Update `multiqc` version from `1.10.1` to `1.11`

### Parameters

- Added `--output_demultiplex_fast5` to output demultiplexed fast5
- Added `--trim_barcodes` in Guppy basecaller to trim the barcodes from output fastq
- Added `--call_variants` to detect DNA variants
- Added `--split_mnps` to split multi-nucleotide polymorphisms into single nucleotide polymorphisms when using medaka
- Added `--phase_vcf` to output a phased vcf when using medaka
- Added `--deepvariant_gpu` to use gpu with docker pepper_margin_deepvariant
- Added `--skip_vc` to skip `variant_calling`
- Added `--skip_sv` to skip `structural_variant_calling`
- Added `--variant_caller` to specify variant caller.
- Added `--structural_variant_caller` to specify structural variant caller
- Added `--skip_modification_analysis` to skip RNA modification detection
- Added `--skip_xpore` to skip `xpore`
- Added `--skip_m6anet` to skip `m6anet`
- Added `--skip_fusion_analysis` to skip RNA fusion detection
- Added `--jaffal_ref_dir` to indicate the reference directory path required by `JAFFAL`

### Software dependencies

| Dependency                      | Old version | New version |
| ------------------------------- | ----------- | ----------- |
| `bioconductor-bambu`            | 1.0.2       | 2.0.0       |
| `bioconductor-bsgenome`         | 1.58.0      | 1.62.0      |
| `cutesv`                        |             | 1.0.12      |
| `deepvariant`                   |             | 1.0.3       |
| `jaffa`                         |             | 2.0         |
| `m6anet`                        |             | 1.0         |
| `medaka`                        |             | 1.4.4       |
| `multiqc`                       | 1.10.1      | 1.11        |
| `ont_fast5_api`                 |             | 4.0.0       |
| `pepper_margin_deepvariant`     |             | 0.8         |
| `pepper_margin_deepvariant_gpu` |             | 0.8         |
| `samtools`                      | 1.14        | 1.15        |
| `sniffles`                      |             | 1.0.12      |
| `xpore`                         |             | 2.1         |

### Bug fix

- The `GET_TEST_DATA` process now uses checks for any file in the path.

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [2.0.1] - 2021-11-29

### Bug fix

- The `UCSC_BEDGRAPHTOBIGWIG` process now uses the `ucsc-bedgraphtobigwig` container
- The full-size and minimal AWS tests have successfully finished after changing to the `ucsc-bedgraphtobigwig` container

## [2.0.0] - 2021-11-26

### Major enhancements

- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
- Software containers are now obtained from [Biocontainers](https://biocontainers.pro/#/registry)
- Update pipeline template to nf-core/tools `2.1`
- [#77](https://github.com/nf-core/nanoseq/issues/77) - Skipped alignment steps
- [#97](https://github.com/nf-core/nanoseq/issues/97) - Add optional DNA cleaning option

### Parameters

- Added `--run_nanolyse` to run NanoLyse for DNA cleaning of FastQ files
- Added `--nanolyse_fasta` to provide a fasta file for nanolyse to filter against

### Software dependencies

| Dependency           | Old version | New version |
| -------------------- | ----------- | ----------- |
| `bioconductor-bambu` | 1.0.0       | 1.0.2       |
| `nanolyse`           |             | 1.2.0       |
| `r-base`             | 4.0.3       | 4.0.2       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [1.1.0] - 2020-11-06

### Major enhancements

- Transcript reconstruction and quantification ([`bambu`](https://bioconductor.org/packages/release/bioc/html/bambu.html) or [`StringTie2`](https://ccb.jhu.edu/software/stringtie/) and [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
- Differential expression analysis at the gene-level ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and transcript-level ([`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html))
- Ability to provide BAM input to the pipeline
- Change samplesheet format to be more flexible to BAM input files
- Add pycoQC and featureCounts output to MultiQC report
- Add AWS full-sized test data
- Add parameter JSON schema for pipeline
- Add citations file
- Update pipeline template to nf-core/tools `1.11`
- Collapsible sections for output files in `docs/output.md`
- Replace `set` with `tuple` and `file` with `path` in `input` section of all processes
- Capitalise process names
- Added `--gpus all` to Docker `runOptions` when using GPU as mentioned [here](https://github.com/docker/compose/issues/6691#issuecomment-514429646)
- Cannot invoke method `containsKey()` on null object when `--igenomes_ignore` is set [#76](https://github.com/nf-core/nanoseq/issues/76)

### Parameters

- Added `--barcode_both_ends` requires barcode on both ends for Guppy basecaller
- Added `--quantification_method` to specify the transcript quantification method to use
- Added `--skip_quantification` to skip transcript quantification and differential analysis
- Added `--skip_differential_analysis` to skip differential analysis with DESeq2 and DEXSeq
- Added `--publish_dir_mode` to customise method of publishing results to output directory [nf-core/tools#585](https://github.com/nf-core/tools/issues/585)

### Software dependencies

| Dependency              | Old version | New version |
| ----------------------- | ----------- | ----------- |
| `Guppy`                 | 3.4.4       | 4.0.14      |
| `markdown`              | 3.1.1       | 3.3.3       |
| `multiqc`               | 1.8         | 1.9         |
| `nanoplot`              | 1.28.4      | 1.32.1      |
| `pygments`              | 2.5.2       | 2.7.2       |
| `pymdown-extensions`    | 6.0         | 8.0.1       |
| `python`                | 3.7.3       | 3.8.6       |
| `samtools`              | 1.9         | 1.11        |
| `ucsc-bedgraphtobigwig` | 357         | 377         |
| `ucsc-bedtobigbed`      | 357         | 377         |
| `bioconductor-bambu`    | -           | 1.0.0       |
| `bioconductor-bsgenome` | -           | 1.58.0      |
| `bioconductor-deseq2`   | -           | 1.30.0      |
| `bioconductor-dexseq`   | -           | 1.36.0      |
| `bioconductor-drimseq`  | -           | 1.18.0      |
| `bioconductor-stager`   | -           | 1.12.0      |
| `r-base`                | -           | 4.0.3       |
| `seaborn`               | -           | 0.10.1      |
| `stringtie`             | -           | 2.1.4       |
| `subread`               | -           | 2.0.1       |
| `psutil`                | -           | -           |

> **NB:** Dependency has been **updated** if both old and new version information is present.
> **NB:** Dependency has been **added** if just the new version information is present.
> **NB:** Dependency has been **removed** if version information isn't present.

## [1.0.0] - 2020-03-05

Initial release of nf-core/nanoseq, created with the [nf-core](http://nf-co.re/) template.
