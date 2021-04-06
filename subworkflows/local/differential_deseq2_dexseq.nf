/*
 * Differential Expression Analysis with DESeq2 and DEXSeq
 */

params.deseq2_options   = [:]
params.dexseq_options   = [:]

include { DESEQ2      } from '../../modules/local/deseq2'       addParams( options: params.deseq2_options )
include { DEXSEQ      } from '../../modules/local/dexseq'       addParams( options: params.dexseq_options )


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
    deseq2_version = DESEQ2.out.deseq2_version
    r_version = DESEQ2.out.r_version
  
    /*
     * DEXseq differential expression of transcripts
     */
    DEXSEQ ( ch_transcript_counts )
    ch_dexseq_txt  = DEXSEQ.out.dexseq_txt
    dexseq_version = DEXSEQ.out.dexseq_version
    drimseq_version = DEXSEQ.out.drimseq_version
    stager_version = DEXSEQ.out.stager_version

    emit:
    ch_deseq2_txt
    ch_dexseq_txt
    deseq2_version
    dexseq_version
    drimseq_version
    stager_version
    r_version
}
