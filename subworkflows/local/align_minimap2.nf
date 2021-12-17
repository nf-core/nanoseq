/*
 * Alignment with MINIMAP2
 */

include { MINIMAP2_INDEX          } from '../../modules/local/minimap2_index'
include { MINIMAP2_ALIGN          } from '../../modules/local/minimap2_align'

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
    minimap2_version = MINIMAP2_INDEX.out.versions

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

    emit:
    ch_index
    minimap2_version
    ch_align_sam
}
