/*
 * Transcript Discovery and Quantification with StringTie2 and FeatureCounts
 */

include { STRINGTIE_STRINGTIE } from '../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE } from '../../modules/nf-core/stringtie/merge/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE } from '../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPT } from '../../modules/nf-core/subread/featurecounts/main'

workflow QUANTIFY_STRINGTIE_FEATURECOUNTS {
    take:
    ch_sorted_bam

    main:

    /*
     * Novel isoform detection with StringTie
     */
    ch_annotation_gtf = Channel.from(file(params.gtf))
    ch_sorted_bam
        .combine(ch_annotation_gtf)
        .map { it -> it[2] }
        .set { ch_annotation_gtf }

    STRINGTIE_STRINGTIE ( ch_sorted_bam, ch_annotation_gtf )
    ch_stringtie_gtf   = STRINGTIE_STRINGTIE.out.transcript_gtf
    stringtie2_version = STRINGTIE_STRINGTIE.out.versions

    /*
     * Merge isoforms across samples called by StringTie
     */
    STRINGTIE_MERGE ( ch_stringtie_gtf.collect{it[1]}, ch_annotation_gtf.unique() )
    ch_stringtie_merged_gtf = STRINGTIE_MERGE.out.gtf

    /*
     * Gene and transcript quantification with featureCounts
     */
    ch_stringtie_merged_gtf
        .combine( [[id:'gene']] )
        .combine( [ch_sorted_bam.collect{it[1]}])
        .map {it -> [it[1], it[2].value, it[0]]}
        .set { ch_featurecounts_gene_input }
    ch_stringtie_merged_gtf
        .combine( [[id:'transcript']] )
        .combine( [ch_sorted_bam.collect{it[1]}])
        .map {it -> [it[1], it[2].value, it[0]]}
        .set { ch_featurecounts_transcript_input }

    SUBREAD_FEATURECOUNTS_GENE ( ch_featurecounts_gene_input )
    SUBREAD_FEATURECOUNTS_TRANSCRIPT ( ch_featurecounts_transcript_input )
    ch_gene_counts                   = SUBREAD_FEATURECOUNTS_GENE.out.counts.map{it -> it[1]}
    ch_transcript_counts             = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.counts.map{it -> it[1]}
    featurecounts_gene_multiqc       = SUBREAD_FEATURECOUNTS_GENE.out.summary
    featurecounts_transcript_multiqc = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.summary
    featurecounts_version            = SUBREAD_FEATURECOUNTS_GENE.out.versions

    emit:
    ch_stringtie_gtf
    ch_stringtie_merged_gtf
    stringtie2_version
    ch_gene_counts
    ch_transcript_counts
    featurecounts_gene_multiqc
    featurecounts_transcript_multiqc
    featurecounts_version
}
