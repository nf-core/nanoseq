/*
 * Differential Expression Analysis with DESeq2 and DEXSeq
 */

params.deseq2_options   = [:]
params.dexseq_options   = [:]

include { DESEQ2      } from '../process/deseq2'       addParams( options: params.deseq2_options )
// include { DEXSEQ      } from '../process/dexseq'       addParams( options: params.dexseq_options )


workflow DIFFERENTIAL_DESEQ2_DEXSEQ {
    take:
    ch_gene_counts
    ch_transcript_counts

    main:
    /*
     * DESeq2 differential expression of genes
     */
    DESEQ2 ( ch_gene_counts )
    ch_deseq2_txt = DESEQ2.out.deseq2_txt
  
    /*
     * DEXseq differential expression of transcripts
     */
 //   DEXSEQ ( ch_transcript_counts )
    
    emit:
    ch_deseq2_txt
//    ch_dexseq_txt
}
