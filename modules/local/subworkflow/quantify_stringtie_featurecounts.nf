/*
 * Transcript Discovery and Quantification with StringTie2 and FeatureCounts
 */

params.stringtie2_options      = [:]
params.featurecounts_options   = [:]

include { STRINGTIE2            } from '../process/stringtie2'             addParams( options: params.stringtie2_options   )
include { STRINGTIE2_MERGE      } from '../process/stringtie2_merge'       addParams( options: params.stringtie2_options    )
include { SUBREAD_FEATURECOUNTS } from '../process/subread_featurecounts'  addParams( options: params.featurecounts_options )

workflow QUANTIFY_STRINGTIE_FEATURECOUNTS {
    take:
    ch_sample
    ch_sortbam

    main:

    ch_sample
        .map  { it -> [ it[0], it[2], it[3] ] }
        .join ( ch_sortbam )
        .set  { ch_sample }

    /*
     * Novel isoform detection with StringTie
     */
    STRINGTIE2 ( ch_sample )
    ch_stringtie_gtf   = STRINGTIE2.out.stringtie_gtf
    stringtie2_version = STRINGTIE2.out.version

    ch_sample
        .map { it -> [ it[2] ] }
        .unique()
        .set { ch_sample_gtf }

    /*
     * Merge isoforms across samples caleed by StringTie
     */
    STRINGTIE2_MERGE ( ch_stringtie_gtf, ch_sample_gtf )
    ch_stringtie_merged_gtf = STRINGTIE2_MERGE.out.merged_gtf
    
    /*
     * Gene and transcript quantification with featureCounts
     */
    ch_sample
        .collect { it[-1]    }
        .set     { ch_sample }
    SUBREAD_FEATURECOUNTS ( ch_stringtie_merged_gtf, ch_sample )
    ch_gene_counts                   = SUBREAD_FEATURECOUNTS.out.gene_counts
    ch_transcript_counts             = SUBREAD_FEATURECOUNTS.out.transcript_counts
    featurecounts_gene_multiqc       = SUBREAD_FEATURECOUNTS.out.featurecounts_gene_multiqc
    featurecounts_transcript_multiqc = SUBREAD_FEATURECOUNTS.out.featurecounts_transcript_multiqc
    featurecounts_version            = SUBREAD_FEATURECOUNTS.out.version

    emit:
    ch_stringtie_gtf
    ch_stringtie_merged_gtf
    ch_gene_counts
    ch_transcript_counts
    featurecounts_gene_multiqc
    featurecounts_transcript_multiqc
    stringtie2_version
    featurecounts_version
}
