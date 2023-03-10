/*
 * Alignment with GRAPHMAP2
 */

include { GRAPHMAP2_INDEX         } from '../../modules/local/graphmap2_index'
include { GRAPHMAP2_ALIGN         } from '../../modules/local/graphmap2_align'

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
        .collate(12) // [fasta, fasta sizes, gtf, bed, fasta_index, annotation_string, meta, fastq, fasta, gtf, is_transcript, fasta_gtf_string]
        .map { it -> [ it[6], it[7], it[0], it[1], it[2], it[3], it[10], it[4] ] } // [ sample, fastq, fasta, sizes, gtf, bed, is_transcripts, index ]
        .set { ch_index }

    /*
     * Map reads with GRAPHMAP2
     */

    GRAPHMAP2_ALIGN ( ch_index )
    ch_align_sam = GRAPHMAP2_ALIGN.out.align_sam

    emit:
    ch_index
    graphmap2_version
    ch_align_sam
}
