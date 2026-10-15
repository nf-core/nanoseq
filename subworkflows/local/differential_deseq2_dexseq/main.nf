/*
 * Differential Expression Analysis with DESeq2 and DEXSeq
 */

include { DESEQ2 } from '../../../modules/local/deseq2/main'
include { DEXSEQ } from '../../../modules/local/dexseq/main'

workflow DIFFERENTIAL_DESEQ2_DEXSEQ {
    take:
    ch_gene_counts        // channel: path(counts)
    ch_transcript_counts  // channel: path(counts)
    quantification_method // val: 'bambu' or 'stringtie2'

    main:
    /*
     * DESeq2 differential expression of genes
     */
    DESEQ2 ( ch_gene_counts, quantification_method )

    /*
     * DEXseq differential expression of transcripts
     */
    DEXSEQ ( ch_transcript_counts, quantification_method )

    emit:
    deseq2_txt = DESEQ2.out.deseq2_txt // channel: path(txt)
    dexseq_txt = DEXSEQ.out.dexseq_txt // channel: path(txt)
}
