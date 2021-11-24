/*
 * Alignment with GRAPHMAP2
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { GRAPHMAP2_INDEX         } from '../../modules/local/graphmap2_index'       addParams( options: params.index_options    )
include { GRAPHMAP2_ALIGN         } from '../../modules/local/graphmap2_align'       addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS       } from './bam_sort_samtools'              addParams( options: params.samtools_options )

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
    graphmap2_version = GRAPHMAP2_INDEX.out.versions

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
    BAM_SORT_SAMTOOLS ( ch_align_sam )
    ch_sortbam               = BAM_SORT_SAMTOOLS.out.sortbam
    ch_sortbam_stats_multiqc = BAM_SORT_SAMTOOLS.out.sortbam_stats_multiqc
    samtools_version         = BAM_SORT_SAMTOOLS.out.versions

    emit:
    graphmap2_version
    ch_sortbam
    ch_sortbam_stats_multiqc
    samtools_version
}
