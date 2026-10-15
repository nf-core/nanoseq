/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules
//
include { GET_TEST_DATA                    } from '../modules/local/get_test_data'
include { GET_NANOLYSE_FASTA               } from '../modules/local/get_nanolyse_fasta'
include { FAST5_TO_POD5                    } from '../modules/local/fast5_to_pod5'
include { DORADO_BASECALLER                } from '../modules/local/dorado_basecaller'
include { DORADO_ALIGNER                   } from '../modules/local/dorado_aligner'
include { MODKIT_PILEUP                    } from '../modules/local/modkit_pileup'
include { BAM_RENAME                       } from '../modules/local/bam_rename'
include { BAMBU                            } from '../modules/local/bambu'

//
// MODULE: Installed directly from nf-core/modules
//
include { QCAT                             } from '../modules/nf-core/qcat/main'
include { NANOLYSE                         } from '../modules/nf-core/nanolyse/main'
include { CUSTOM_GETCHROMSIZES             } from '../modules/nf-core/custom/getchromsizes/main'
include { SAMTOOLS_INDEX                   } from '../modules/nf-core/samtools/index/main'
include { MULTIQC                          } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOW: Local subworkflows
//
include { INPUT_CHECK                      } from '../subworkflows/local/input_check'
include { QCFASTQ_NANOPLOT_FASTQC          } from '../subworkflows/local/qcfastq_nanoplot_fastqc'
include { ALIGN_MINIMAP2                   } from '../subworkflows/local/align_minimap2'
include { ALIGN_GRAPHMAP2                  } from '../subworkflows/local/align_graphmap2'
include { BEDTOOLS_UCSC_BIGWIG             } from '../subworkflows/local/bedtools_ucsc_bigwig'
include { BEDTOOLS_UCSC_BIGBED             } from '../subworkflows/local/bedtools_ucsc_bigbed'
include { QUANTIFY_STRINGTIE_FEATURECOUNTS } from '../subworkflows/local/quantify_stringtie_featurecounts'
include { DIFFERENTIAL_DESEQ2_DEXSEQ       } from '../subworkflows/local/differential_deseq2_dexseq'
include { RNA_MODIFICATION_XPORE_M6ANET    } from '../subworkflows/local/rna_modifications_xpore_m6anet'
include { RNA_FUSIONS_JAFFAL               } from '../subworkflows/local/rna_fusions_jaffal'

//
// FUNCTIONS
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_nanoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOSEQ {

    take:
    ch_samplesheet              // channel: path to samplesheet read in from --input
    profile                     // string : profile(s) the pipeline was launched with
    multiqc_config              // string : path to custom MultiQC config
    multiqc_logo                // string : path to custom MultiQC logo
    multiqc_methods_description // string : path to custom MultiQC methods description
    outdir                      // string : output directory

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    def rna_protocol = params.protocol == 'cDNA' || params.protocol == 'directRNA'
    def fasta        = params.fasta ? file(params.fasta, checkIfExists: true) : null
    def gtf          = params.gtf   ? file(params.gtf, checkIfExists: true)   : null

    //
    // Resolve raw data directory given with --input_path
    //
    // Nextflow cannot recursively download directories over HTTPS, so for the test profiles
    // that need raw FAST5/POD5 directories the nf-core/test-datasets repository is cloned.
    //
    def profiles        = profile.tokenize(',')
    def is_test_profile = profiles.any { p -> p.startsWith('test') }
    def nonbarcoded     = profiles.any { p -> p == 'test_bc_nodx' || p.contains('rnamod') }

    def ch_input_path
    if (is_test_profile && (!params.skip_basecalling || !params.skip_modification_analysis)) {
        GET_TEST_DATA ( nonbarcoded ? 'nonbarcoded_multi' : 'barcoded_multi' )
        if (!params.skip_modification_analysis) {
            ch_input_path = GET_TEST_DATA.out.modification_dir
        } else if (params.bedmethyl_out) {
            ch_input_path = GET_TEST_DATA.out.pod5_dir
        } else {
            ch_input_path = GET_TEST_DATA.out.fast5_dir
        }
    } else if (params.input_path) {
        ch_input_path = channel.fromPath(params.input_path, checkIfExists: true)
    } else {
        ch_input_path = channel.value('not_changed')
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // Relative 'input_file' entries are only resolved against the cloned test data
    //
    INPUT_CHECK (
        ch_samplesheet,
        is_test_profile && nonbarcoded ? ch_input_path : channel.value('not_changed')
    )
    def ch_sample = INPUT_CHECK.out.sample

    def ch_fastq = channel.empty()
    def barcode_kit = params.barcode_kit ?: 'Auto'

    if (!params.skip_basecalling) {

        def ch_raw = params.skip_demultiplexing
            ? ch_sample
            : ch_input_path.map { path -> [ [ id: 'undemultiplexed' ], path ] }

        def ch_pod5
        if (params.input_path_file_type == 'fast5') {
            //
            // MODULE: Convert FAST5 to POD5
            //
            FAST5_TO_POD5 ( ch_raw )
            ch_pod5 = FAST5_TO_POD5.out.pod5
        } else {
            ch_pod5 = ch_raw
        }

        //
        // MODULE: Basecalling (and modified base calling) with dorado
        //
        DORADO_BASECALLER (
            ch_pod5,
            params.dorado_device,
            params.dorado_model,
            params.bedmethyl_out ? params.dorado_modification : null
        )

        if (params.bedmethyl_out) {
            //
            // MODULE: Align modified base calls and create bedMethyl files with modkit
            //
            DORADO_ALIGNER ( DORADO_BASECALLER.out.bam, fasta )
            MODKIT_PILEUP ( DORADO_ALIGNER.out.bam.join(DORADO_ALIGNER.out.bai) )
        } else if (!params.skip_demultiplexing) {
            //
            // MODULE: Demultiplexing with qcat
            //
            QCAT ( DORADO_BASECALLER.out.fastq, barcode_kit )
            ch_fastq = assignBarcodedReads(QCAT.out.reads, ch_sample)
        } else {
            ch_fastq = DORADO_BASECALLER.out.fastq
        }

    } else if (!params.skip_demultiplexing) {

        //
        // MODULE: Demultiplexing with qcat
        //
        QCAT (
            ch_input_path.map { path -> [ [ id: 'undemultiplexed' ], path ] },
            barcode_kit
        )
        ch_fastq = assignBarcodedReads(QCAT.out.reads, ch_sample)

    } else if (!params.skip_alignment || !params.skip_fusion_analysis) {
        ch_fastq = ch_sample.filter { _meta, reads -> reads && reads.name.endsWith('.gz') }
    }

    //
    // MODULE: DNA contaminant removal using NanoLyse
    //
    def ch_fastq_to_align = ch_fastq
    if (params.run_nanolyse) {
        def ch_nanolyse_fasta
        if (params.nanolyse_fasta) {
            ch_nanolyse_fasta = channel.value(file(params.nanolyse_fasta, checkIfExists: true))
        } else {
            GET_NANOLYSE_FASTA ()
            ch_nanolyse_fasta = GET_NANOLYSE_FASTA.out.fasta
        }
        NANOLYSE ( ch_fastq, ch_nanolyse_fasta )
        ch_fastq_to_align = NANOLYSE.out.fastq
    }

    //
    // SUBWORKFLOW: Read QC with NanoPlot, ToulligQC and FastQC
    //
    if (!params.skip_qc) {
        QCFASTQ_NANOPLOT_FASTQC (
            ch_fastq,
            params.skip_nanoplot,
            params.skip_toulligqc,
            params.skip_fastqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(QCFASTQ_NANOPLOT_FASTQC.out.multiqc_files)
    }

    //
    // MODULE: Make chromosome sizes and FASTA index
    //
    def ch_fasta_fai = channel.empty()
    def ch_chr_sizes = channel.empty()
    if (!params.bedmethyl_out) {
        def ch_fasta = channel.value([ [ id: 'reference' ], fasta ])
        CUSTOM_GETCHROMSIZES ( ch_fasta )
        ch_chr_sizes = CUSTOM_GETCHROMSIZES.out.sizes
        ch_fasta_fai = ch_fasta
            .combine(CUSTOM_GETCHROMSIZES.out.fai)
            .map { meta, fa, _meta2, fai -> [ meta, fa, fai ] }
            .first()
    }

    def ch_sorted_bam = channel.empty()
    def ch_sorted_bai = channel.empty()
    if (!params.skip_alignment) {

        if (params.aligner == 'minimap2') {
            //
            // SUBWORKFLOW: Align with minimap2
            //
            ALIGN_MINIMAP2 ( ch_fasta_fai, ch_fastq_to_align )
            ch_sorted_bam = ALIGN_MINIMAP2.out.bam
            ch_sorted_bai = ALIGN_MINIMAP2.out.bai
            ch_multiqc_files = ch_multiqc_files.mix(
                ALIGN_MINIMAP2.out.stats.mix(ALIGN_MINIMAP2.out.flagstat, ALIGN_MINIMAP2.out.idxstats).map { _meta, file -> file }
            )
        } else {
            //
            // SUBWORKFLOW: Align with GraphMap2
            //
            ALIGN_GRAPHMAP2 ( ch_fasta_fai, ch_fastq_to_align )
            ch_sorted_bam = ALIGN_GRAPHMAP2.out.bam
            ch_sorted_bai = ALIGN_GRAPHMAP2.out.bai
            ch_multiqc_files = ch_multiqc_files.mix(
                ALIGN_GRAPHMAP2.out.stats.mix(ALIGN_GRAPHMAP2.out.flagstat, ALIGN_GRAPHMAP2.out.idxstats).map { _meta, file -> file }
            )
        }

        //
        // SUBWORKFLOW: Create bigWig coverage tracks
        //
        if (!params.skip_bigwig) {
            BEDTOOLS_UCSC_BIGWIG ( ch_sorted_bam, ch_chr_sizes )
        }

        //
        // SUBWORKFLOW: Create bigBed coverage tracks
        //
        if (!params.skip_bigbed) {
            BEDTOOLS_UCSC_BIGBED ( ch_sorted_bam, ch_chr_sizes )
        }

    } else {

        //
        // MODULE: Rename pre-aligned BAM files provided in the samplesheet and index them
        //
        BAM_RENAME ( ch_sample.filter { _meta, reads -> reads && reads.name.endsWith('.bam') } )
        ch_sorted_bam = BAM_RENAME.out.bam
        SAMTOOLS_INDEX ( ch_sorted_bam )
        ch_sorted_bai = SAMTOOLS_INDEX.out.index
    }

    //
    // Transcript quantification and differential analysis
    //
    if (!params.skip_quantification && rna_protocol) {

        def ch_gene_counts
        def ch_transcript_counts
        if (params.quantification_method == 'bambu') {
            //
            // MODULE: Quantification and novel isoform detection with bambu
            //
            BAMBU (
                channel.value([ fasta, gtf ]),
                ch_sorted_bam.map { _meta, bam -> bam }.collect()
            )
            ch_gene_counts       = BAMBU.out.gene_counts
            ch_transcript_counts = BAMBU.out.transcript_counts
        } else {
            //
            // SUBWORKFLOW: Novel isoform detection with StringTie2 and quantification with featureCounts
            //
            QUANTIFY_STRINGTIE_FEATURECOUNTS ( ch_sorted_bam, channel.value(gtf) )
            ch_gene_counts       = QUANTIFY_STRINGTIE_FEATURECOUNTS.out.gene_counts
            ch_transcript_counts = QUANTIFY_STRINGTIE_FEATURECOUNTS.out.transcript_counts
            ch_multiqc_files     = ch_multiqc_files.mix(
                QUANTIFY_STRINGTIE_FEATURECOUNTS.out.featurecounts_gene_summary,
                QUANTIFY_STRINGTIE_FEATURECOUNTS.out.featurecounts_transcript_summary
            )
        }

        //
        // SUBWORKFLOW: Differential gene and transcript analysis with DESeq2 and DEXSeq
        //
        if (!params.skip_differential_analysis) {
            DIFFERENTIAL_DESEQ2_DEXSEQ ( ch_gene_counts, ch_transcript_counts, params.quantification_method )
        }
    }

    //
    // SUBWORKFLOW: RNA modification detection with xPore and m6anet
    //
    if (!params.skip_modification_analysis && params.protocol == 'directRNA') {
        def ch_nanopolish_input = ch_sorted_bam
            .join(ch_sorted_bai)
            .join(ch_sample)
            .map { meta, bam, bai, reads -> [ meta, fasta, gtf, reads, bam, bai ] }

        RNA_MODIFICATION_XPORE_M6ANET ( ch_nanopolish_input, params.skip_xpore, params.skip_m6anet )
    }

    //
    // SUBWORKFLOW: RNA fusion detection with JAFFAL
    //
    if (!params.skip_fusion_analysis && rna_protocol) {
        RNA_FUSIONS_JAFFAL ( ch_fastq, params.jaffal_ref_dir )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'nanoseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'nanoseq'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Assign demultiplexed qcat reads (e.g. barcode01.fastq.gz) to samples using the samplesheet barcode
//
def assignBarcodedReads(ch_qcat_reads, ch_sample) {
    return ch_qcat_reads
        .flatMap { _meta, reads -> (reads instanceof List) ? reads : [ reads ] }
        .map { reads -> [ reads.name.replaceAll(/\.f(ast)?q\.gz$/, ''), reads ] }
        .join( ch_sample.map { meta, _reads -> [ meta.barcode, meta ] } )
        .map { _barcode, reads, meta -> [ meta, reads ] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
