/*
 * Convert BAM to BigBed
 */

params.bigbed_options   = [:]

include { BEDTOOLS_BAMBED     } from '../../modules/local/bedtools_bamtobed'  addParams( options: params.bigbed_options )
include { UCSC_BED12TOBIGBED  } from '../../modules/local/ucsc_bed12tobigbed' addParams( options: params.bigbed_options )

workflow BEDTOOLS_UCSC_BIGBED {
    take:
    ch_sortbam
    
    main:
    /*
     * Convert BAM to BED12
     */
    BEDTOOLS_BAMBED ( ch_sortbam )
    ch_bed12         = BEDTOOLS_BAMBED.out.bed12
    bedtools_version = BEDTOOLS_BAMBED.out.version

    /*
     * Convert BED12 to BigBED
     */
    UCSC_BED12TOBIGBED ( ch_bed12 )
    ch_bigbed = UCSC_BED12TOBIGBED.out.bigbed


    emit:
    bedtools_version
    ch_bigbed
}
