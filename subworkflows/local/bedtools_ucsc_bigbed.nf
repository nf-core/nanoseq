/*
 * Convert BAM to BigBed
 */

include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main'
include { UCSC_BEDTOBIGBED } from '../../modules/nf-core/ucsc/bedtobigbed/main'

workflow BEDTOOLS_UCSC_BIGBED {
    take:
    ch_sorted_bam
    ch_chr_sizes

    main:
    /*
     * Convert BAM to BED12
     */
    BEDTOOLS_BAMTOBED ( ch_sorted_bam )
    ch_bed         = BEDTOOLS_BAMTOBED.out.bed
    bedtools_version = BEDTOOLS_BAMTOBED.out.versions

    /*
     * Convert BED12 to BigBED
     */
    ch_bed
        .combine(ch_chr_sizes)
        .map { it -> it[3] }
        .set { ch_sizes }
    UCSC_BEDTOBIGBED ( ch_bed, ch_sizes, [] )
    ch_bigbed = UCSC_BEDTOBIGBED.out.bigbed
    bed12tobigbed_version = UCSC_BEDTOBIGBED.out.versions

    emit:
    ch_bigbed
    bedtools_version
    bed12tobigbed_version
}
