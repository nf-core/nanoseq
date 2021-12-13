/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

params.options = [:]

include { SAMTOOLS_VIEW_BAM   } from '../../modules/local/samtools_view_bam'              addParams( options: params.options )
include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index'            addParams( options: params.options )
include { BAM_STATS_SAMTOOLS  } from '../../subworkflows/nf-core/bam_stats_samtools'      addParams( options: params.options )

workflow BAM_SORT_INDEX_SAMTOOLS {
    take:
    ch_sam // channel: [ val(meta), [ bam ] ]

    main:
    SAMTOOLS_VIEW_BAM  ( ch_sam )

    SAMTOOLS_SORT_INDEX( SAMTOOLS_VIEW_BAM.out.bam )

    ch_sam
        .join( SAMTOOLS_SORT_INDEX.out.bam_bai )
        .map { it -> [ it[0], it[1], it[2], it[4], it[5] ] }
        .set { sortbam }

    BAM_STATS_SAMTOOLS ( SAMTOOLS_SORT_INDEX.out.bam_bai )
    BAM_STATS_SAMTOOLS.out.stats
        .join ( BAM_STATS_SAMTOOLS.out.idxstats )
        .join ( BAM_STATS_SAMTOOLS.out.flagstat )
        .map  { it -> [ it[1], it[2], it[3] ] }
        .set  { sortbam_stats_multiqc }
    versions = BAM_STATS_SAMTOOLS.out.versions

    emit:
    sortbam
    sortbam_stats_multiqc
    versions
}
