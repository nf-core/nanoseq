/*
 * Convert BAM to BigBed
 */

include { BEDTOOLS_BAMBED     } from '../../modules/local/bedtools_bamtobed'
include { UCSC_BED12TOBIGBED  } from '../../modules/local/ucsc_bed12tobigbed'

workflow BEDTOOLS_UCSC_BIGBED {
    take:
    ch_sortbam

    main:
    /*
     * Convert BAM to BED12
     */
    BEDTOOLS_BAMBED ( ch_sortbam )
    ch_bed12         = BEDTOOLS_BAMBED.out.bed12
    bedtools_version = BEDTOOLS_BAMBED.out.versions

    /*
     * Convert BED12 to BigBED
     */
    UCSC_BED12TOBIGBED ( ch_bed12 )
    ch_bigbed = UCSC_BED12TOBIGBED.out.bigbed
    bed12tobigbed_version = UCSC_BED12TOBIGBED.out.versions

    emit:
    bedtools_version
    ch_bigbed
    bed12tobigbed_version
}
