/*
 * Transcript discovery and quantification with StringTie2 and featureCounts
 */

include { STRINGTIE_STRINGTIE                                        } from '../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE                                            } from '../../modules/nf-core/stringtie/merge/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE       } from '../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPT } from '../../modules/nf-core/subread/featurecounts/main'

workflow QUANTIFY_STRINGTIE_FEATURECOUNTS {
    take:
    ch_bam // channel: [ val(meta), path(bam) ]
    ch_gtf // channel (value): path(gtf)

    main:
    /*
     * Novel isoform detection with StringTie (long-read BAMs only, guided by the annotation)
     */
    STRINGTIE_STRINGTIE(
        ch_bam.map { meta, bam -> [meta, [], bam] },
        'long-reads-assembly',
        ch_gtf,
    )

    /*
     * Merge isoforms across samples called by StringTie
     */
    ch_stringtie_gtfs = STRINGTIE_STRINGTIE.out.transcript_gtf
        .map { _meta, gtf -> gtf }
        .collect()
        .map { gtfs -> [[id: 'stringtie.merged'], gtfs] }
    STRINGTIE_MERGE(
        ch_stringtie_gtfs,
        ch_gtf.map { gtf -> [[id: 'annotation'], gtf] },
    )
    ch_merged_gtf = STRINGTIE_MERGE.out.merged_gtf.map { _meta, gtf -> gtf }

    /*
     * Gene and transcript quantification with featureCounts across all samples
     */
    ch_bams = ch_bam
        .map { _meta, bam -> bam }
        .collect()

    SUBREAD_FEATURECOUNTS_GENE(
        ch_bams
            .map { bams -> [[id: 'gene', single_end: true], bams] }
            .combine(ch_merged_gtf)
    )
    SUBREAD_FEATURECOUNTS_TRANSCRIPT(
        ch_bams
            .map { bams -> [[id: 'transcript', single_end: true], bams] }
            .combine(ch_merged_gtf)
    )

    emit:
    stringtie_gtf                    = STRINGTIE_STRINGTIE.out.transcript_gtf                               // channel: [ val(meta), path(gtf) ]
    merged_gtf                       = ch_merged_gtf                                                        // channel: path(gtf)
    gene_counts                      = SUBREAD_FEATURECOUNTS_GENE.out.counts.map { _meta, counts -> counts }        // channel: path(tsv)
    transcript_counts                = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.counts.map { _meta, counts -> counts }  // channel: path(tsv)
    featurecounts_gene_summary       = SUBREAD_FEATURECOUNTS_GENE.out.summary.map { _meta, summary -> summary }       // channel: path(summary)
    featurecounts_transcript_summary = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.summary.map { _meta, summary -> summary } // channel: path(summary)
}
