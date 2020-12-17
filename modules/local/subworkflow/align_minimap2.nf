/*
 * Alignment with MINIMAP2
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { MINIMAP2_INDEX          } from '../process/minimap2_index'                       addParams( options: params.index_options    )
include { MINIMAP2_ALIGN          } from '../process/minimap2_align'                       addParams( options: params.align_options    )
include { SAMTOOLS_BAM_SORT_STATS } from '../process/samtools_bam_sort_stats'              addParams( options: params.samtools_options )

workflow ALIGN_MINIMAP2 {
    take:
    ch_fasta_index // channel: [ val(meta), [ reads ] ]
    ch_fastq_alignment
    
    main:
    /*
     * Create genome/transcriptome index
     */
    MINIMAP2_INDEX ( ch_fasta_index )
    ch_index = MINIMAP2_INDEX.out.index

    ch_index
        .cross(ch_fastq_alignment) { it -> it[-1] }
        .flatten()
        .collate(13)
        .map { it -> [ it[7], it[8], it[0], it[1], it[2], it[3], it[4], it[5] ] } // [ sample, fastq, fasta, sizes, gtf, bed, is_transcripts, index ]
        .set { ch_index }

    /*
     * Map reads with MINIMAP2
     */
    MINIMAP2_ALIGN ( ch_index )
    ch_align_sam = MINIMAP2_ALIGN.out.align_sam

    /*
     * Convert SAM to BAM, sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    SAMTOOLS_BAM_SORT_STATS ( ch_align_sam )
    ch_sortbam               = SAMTOOLS_BAM_SORT_STATS.out.sortbam
    ch_sortbam_stats_multiqc = SAMTOOLS_BAM_SORT_STATS.out.sortbam_stats_multiqc

    emit:
    ch_sortbam
    ch_sortbam_stats_multiqc
//    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
