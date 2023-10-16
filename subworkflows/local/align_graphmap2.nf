/*
 * Alignment with GRAPHMAP2
 */

include { GRAPHMAP2_INDEX } from '../../modules/nf-core/graphmap2/index/main'
include { GRAPHMAP2_ALIGN } from '../../modules/nf-core/graphmap2/align/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools/main'

workflow ALIGN_GRAPHMAP2 {
    take:
    ch_fasta
    ch_fastq

    main:
    /*
     * Create genome/transcriptome index
     */
    ch_fasta
        .map { it -> it[1] }
        .set { ch_fasta_graphmapindex }
    GRAPHMAP2_INDEX ( ch_fasta_graphmapindex )
    ch_graphmap_index = GRAPHMAP2_INDEX.out.index
    graphmap2_version = GRAPHMAP2_INDEX.out.versions

    /*
     * Map reads with GRAPHMAP2
     */
    ch_fastq
        .map { it -> [ it[0], it[1] ] }
        .set { ch_alignment_input }
    ch_alignment_input
        .combine (ch_fasta_graphmapindex)
        .map { it -> it[2] }
        .set { ch_reference }
    ch_alignment_input
        .combine (ch_graphmap_index)
        .map { it -> it[2] }
        .set { ch_reference_index }

    GRAPHMAP2_ALIGN ( ch_alignment_input, ch_reference, ch_reference_index )
    ch_alignment_input
        .map { it -> it[1] }
        .set { ch_just_fastq }
    GRAPHMAP2_ALIGN.out.sam
        .combine (ch_graphmap_index)
        .combine (ch_just_fastq)
        .map { it -> [ it[0], it[1], it[2], it[3] ] }
        .set { ch_merged_input }
    ch_merged_input
        .map { it -> it[2] }
        .set { ch_notneeded_qname }
    ch_merged_input
        .map { it -> [ it[0], it[1], it[3] ]}
        .set { ch_samtools_input }
    ch_samtools_input.view()
    ch_notneeded_qname.view()
    SAMTOOLS_VIEW ( ch_samtools_input, ch_fasta, ch_notneeded_qname )
    SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
    ch_sorted_bam = SAMTOOLS_SORT.out.bam
    SAMTOOLS_INDEX ( ch_sorted_bam )
    ch_sorted_bai = SAMTOOLS_INDEX.out.bai
    samtools_version=SAMTOOLS_INDEX.out.versions

    ch_sorted_bam
        .join(ch_sorted_bai, by: 0)
        .set { ch_bam_bai }
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_stats = BAM_STATS_SAMTOOLS.out.stats
    ch_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    ch_idxstats = BAM_STATS_SAMTOOLS.out.idxstats

    emit:
    ch_graphmap_index
    graphmap2_version
    ch_sorted_bam
    ch_sorted_bai
    samtools_version
    ch_stats
    ch_flagstat
    ch_idxstats
}
