/*
 * Alignment with MINIMAP2
 */

include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_OTHER } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_VARIANT } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow ALIGN_MINIMAP2 {
    take:
    ch_fasta
    ch_fastq

    main:
    /*
     * Create genome/transcriptome index
     */
    MINIMAP2_INDEX ( ch_fasta )
    ch_minimap_index   = MINIMAP2_INDEX.out.index
    minimap2_version = MINIMAP2_INDEX.out.versions

    /*
     * Map reads with MINIMAP2
     */
    ch_fastq
        .map { it -> [ it[0], it[1] ] }
        .set { ch_alignment_input }
    ch_alignment_input
        .combine ( ch_minimap_index )
        .map { it -> it[3] }
        .set { ch_reference }
    bam_format = true
    cigar_paf_format = false
    cigar_bam = false
    if (params.call_variants) {
        MINIMAP2_ALIGN_VARIANT ( ch_alignment_input, ch_reference, bam_format, cigar_paf_format, cigar_bam )
        ch_sorted_bam = MINIMAP2_ALIGN_VARIANT.out.bam
    } else {
        MINIMAP2_ALIGN_OTHER ( ch_alignment_input, ch_reference, bam_format, cigar_paf_format, cigar_bam )
        ch_sorted_bam = MINIMAP2_ALIGN_OTHER.out.bam
    }

    SAMTOOLS_INDEX ( ch_sorted_bam )
    ch_sorted_bai = SAMTOOLS_INDEX.out.bai
    samtools_version = SAMTOOLS_INDEX.out.versions

    ch_sorted_bam
        .join(ch_sorted_bai, by: 0)
        .set { ch_bam_bai }
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_stats = BAM_STATS_SAMTOOLS.out.stats
    ch_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    ch_idxstats = BAM_STATS_SAMTOOLS.out.idxstats

    emit:
    ch_minimap_index
    minimap2_version
    ch_sorted_bam
    ch_sorted_bai
    samtools_version
    ch_stats
    ch_flagstat
    ch_idxstats
}

