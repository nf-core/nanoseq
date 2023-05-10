/*
 * Convert BAM to BigWig
 */

include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDTOOLS_UCSC_BIGWIG {
    take:
    ch_sorted_bam
    ch_chr_sizes

    main:
    /*
     * Convert BAM to BEDGraph
     */
    ch_sorted_bam
        .combine([1])
        .set { ch_genomecov_input }
    ch_genomecov_input
        .combine(ch_chr_sizes)
        .map { it -> it[4] }
        .set { ch_sizes }
    extension = 'bedGraph'

    BEDTOOLS_GENOMECOV ( ch_genomecov_input, ch_sizes, extension )
    ch_bedgraph      = BEDTOOLS_GENOMECOV.out.genomecov
    bedtools_version = BEDTOOLS_GENOMECOV.out.versions

    /*
     * Convert BEDGraph to BigWig
     */
    UCSC_BEDGRAPHTOBIGWIG ( ch_bedgraph, ch_sizes )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bedgraphtobigwig_version = UCSC_BEDGRAPHTOBIGWIG.out.versions

    emit:
    ch_bigwig
    bedtools_version
    bedgraphtobigwig_version
}
