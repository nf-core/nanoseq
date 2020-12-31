/*
 * Alignment with GRAPHMAP2
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { GRAPHMAP2_INDEX         } from '../process/graphmap2_index'                     addParams( options: params.index_options    )
include { GRAPHMAP2_ALIGN         } from '../process/graphmap2_align'                     addParams( options: params.align_options    )
include { SAMTOOLS_BAM_SORT_STATS } from '../process/samtools_bam_sort_stats'             addParams( options: params.samtools_options )

workflow ALIGN_GRAPHMAP2 {
    take:
    ch_fasta_index // channel: [ val(meta), [ reads ] ]
    ch_fastq
    
    main:
    /*
     * Create genome/transcriptome index
     */
    GRAPHMAP2_INDEX ( ch_fasta_index )
    ch_index          = GRAPHMAP2_INDEX.out.index
    graphmap2_version = GRAPHMAP2_INDEX.out.version 

    ch_index
        .cross(ch_fastq) { it -> it[-1] }
        .flatten()
        .collate(13)
        .map { it -> [ it[6], it[7], it[0], it[1], it[2], it[3], it[10], it[4] ] } // [ sample, fastq, fasta, sizes, gtf, bed, is_transcripts, index ]
        .set { ch_index }

    /*
     * Map reads with GRAPHMAP2
     */
    GRAPHMAP2_ALIGN ( ch_index )
    ch_align_sam = GRAPHMAP2_ALIGN.out.align_sam

    /*
     * Convert SAM to BAM, sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    SAMTOOLS_BAM_SORT_STATS ( ch_align_sam )
    ch_sortbam               = SAMTOOLS_BAM_SORT_STATS.out.sortbam
    ch_sortbam_stats_multiqc = SAMTOOLS_BAM_SORT_STATS.out.sortbam_stats_multiqc
    samtools_version         = SAMTOOLS_BAM_SORT_STATS.out.version

    emit:
    graphmap2_version
    ch_sortbam
    ch_sortbam_stats_multiqc
    samtools_version
}
