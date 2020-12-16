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
    ch_sample_annotation
    ch_sortbam

    main:
    ch_sample_annotation
            .map  { it -> [ it[0], it[2], it[3] ] }
            .join ( ch_sortbam )
            .into { 
                ch_sample_annotation
                ch_sample_gtf
                ch_sample_featurecounts
            }

    /*
     * Novel isoform detection with StringTie
     */
    ch_stringtie_gtf        = STRINGTIE2 ( ch_sample_annotation )

    ch_sample_gtf
            .map { it -> [ it[2]] }
            .unique()
    
    /*
     * Merge isoforms across samples caleed by StringTie
     */
    ch_stringtie_merged_gtf = STRINGTIE2_MERGE ( ch_stringtie_gtf, ch_sample_gtf )

    /*
     * Gene and transcript quantification with featureCounts
     */
    SUBREAD_FEATURECOUNTS ( ch_stringtie_merged_gtf, ch_sample_featurecounts )
    
    emit:
    ch_gene_counts                      = "counts_gene.txt"
    ch_transcript_counts                = "counts_transcript.txt"
    ch_featurecounts_gene_multiqc       = "counts_gene.txt.summary"
    ch_featurecounts_transcript_multiqc = "counts_transcript.txt.summary"
}
