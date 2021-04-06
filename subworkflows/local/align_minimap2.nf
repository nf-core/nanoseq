/*
 * Alignment with MINIMAP2
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { MINIMAP2_INDEX          } from '../../modules/local/minimap2_index'        addParams( options: params.index_options    )
include { MINIMAP2_ALIGN          } from '../../modules/local/minimap2_align'        addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS       } from './bam_sort_samtools'              addParams( options: params.samtools_options )

workflow ALIGN_MINIMAP2 {
    take:
    ch_fasta_index // channel: [ val(meta), [ reads ] ]
    ch_fastq
    
    main:
    /*
     * Create genome/transcriptome index
     */
    MINIMAP2_INDEX ( ch_fasta_index )
    ch_index         = MINIMAP2_INDEX.out.index
    minimap2_version = MINIMAP2_INDEX.out.version

    ch_index
        .cross(ch_fastq) { it -> it[-1] }
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
    BAM_SORT_SAMTOOLS ( ch_align_sam )
    ch_sortbam               = BAM_SORT_SAMTOOLS.out.sortbam
    ch_sortbam_stats_multiqc = BAM_SORT_SAMTOOLS.out.sortbam_stats_multiqc
    samtools_version         = BAM_SORT_SAMTOOLS.out.version

    emit:
    minimap2_version

    ch_sortbam
    ch_sortbam_stats_multiqc
    samtools_version
}
