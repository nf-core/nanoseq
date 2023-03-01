/*
 * Alignment with GRAPHMAP2
 */

include { GRAPHMAP2_INDEX as GRAPHMAP2_INDEX_VARIANT } from '../../modules/nf-core/graphmap2/index/main'
include { GRAPHMAP2_ALIGN as GRAPHMAP2_ALIGN_VARIANT } from '../../modules/nf-core/graphmap2/align/main'
//include { GRAPHMAP2_INDEX as GRAPHMAP2_INDEX_SPLICE } from '../../modules/nf-core/graphmap2/index/main'
//include { GRAPHMAP2_ALIGN as GRAPHMAP2_ALIGN_SPLICE } from '../../modules/nf-core/graphmap2/align/main'

workflow ALIGN_GRAPHMAP2 {
    take:
    ch_fastq // channel: [ val(meta), path(reads)  ]
    ch_fasta // channel: [ path fasta ]
    ch_fai // channel: [ path fai ]

    main:

    ch_versions = Channel.empty()

    if (params.protocol == "DNA") {

        /*
        * Create genome/transcriptome index for variant calling with GRAPHMAP2
        */
        GRAPHMAP2_INDEX_VARIANT (ch_fasta)
        ch_index          = GRAPHMAP2_INDEX_VARIANT.out.index
        ch_versions = ch_versions.mix(GRAPHMAP2_INDEX_VARIANT.out.versions)

        /*
        * Map reads with GRAPHMAP2
        */

        GRAPHMAP2_ALIGN_VARIANT ( ch_index )
        ch_align_sam = GRAPHMAP2_ALIGN.out.align_sam

    }

    emit:
    ch_index
    graphmap2_version
    ch_align_sam
}
