/*
 * Alignment with GRAPHMAP2
 */

params.index_options    = [:]
params.align_options    = [:]
params.samtools_options = [:]

include { GRAPHMAP2_INDEX      } from '../process/graphmap2_index'                     addParams( options: params.index_options    )
include { GRAPHMAP2_ALIGN      } from '../process/graphmap2_align'                     addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS   } from '../../nf-core/subworkflow/bam_sort_samtools'     addParams( options: params.samtools_options )

workflow ALIGN_GRAPHMAP2 {
    take:
    ch_fasta_index // channel: [ val(meta), [ reads ] ]
    ch_fastq_alignment
    
    main:
    /*
     * Create genome/transcriptome index
     */
    ch_index = GRAPHMAP2_INDEX ( ch_fasta_index )

    ch_index
        .cross(ch_fastq_alignment) { it -> it[-1] }
        .flatten()
        .collate(13)
        .map { it -> [ it[7], it[8], it[0], it[1], it[2], it[3], it[4], it[5] ] } // [ sample, fastq, fasta, sizes, gtf, bed, is_transcripts, index ]
        .set { ch_index }

    /*
     * Map reads with GRAPHMAP2
     */
    ch_align_sam = GRAPHMAP2_ALIGN ( ch_index )

    /*
     * Sort, index BAM file and run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( ch_align_sam )

    emit:
    bam              = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}
