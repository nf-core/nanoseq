/*
 * Differential Expression Analysis with DESeq2 and DEXSeq
 */

include { DESEQ2      } from '../../modules/local/deseq2'
include { DEXSEQ      } from '../../modules/local/dexseq'

workflow DIFFERENTIAL_DESEQ2_DEXSEQ {
    take:
    ch_gene_counts
    ch_transcript_counts

    main:
    /*
     * DESeq2 differential expression of genes
     */
    DESEQ2 ( ch_gene_counts )
    ch_deseq2_txt  = DESEQ2.out.deseq2_txt
    deseq2_version = DESEQ2.out.versions

    /*
     * DEXseq differential expression of transcripts
     */
    DEXSEQ ( ch_transcript_counts )
    ch_dexseq_txt  = DEXSEQ.out.dexseq_txt
    dexseq_version = DEXSEQ.out.versions

    emit:
    ch_deseq2_txt
    ch_dexseq_txt
    deseq2_version
    dexseq_version
}
