/*
 * Convert BAM to BigBed
 */

include { BEDTOOLS_BAMTOBED } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_SORT     } from '../../../modules/nf-core/bedtools/sort/main'
include { UCSC_BEDTOBIGBED  } from '../../../modules/nf-core/ucsc/bedtobigbed/main'

workflow BEDTOOLS_UCSC_BIGBED {
    take:
    ch_bam   // channel: [ val(meta), path(bam) ]
    ch_sizes // channel (value): [ val(meta), path(sizes) ]

    main:
    /*
     * Convert BAM to BED12
     */
    BEDTOOLS_BAMTOBED(ch_bam)

    /*
     * Sort BED12 by chromosome and start position as required by bedToBigBed
     */
    BEDTOOLS_SORT(BEDTOOLS_BAMTOBED.out.bed, [], 'bed')

    /*
     * Convert BED12 to BigBed
     */
    UCSC_BEDTOBIGBED(
        BEDTOOLS_SORT.out.sorted,
        ch_sizes.map { _meta, sizes -> sizes },
        [], // autosql
    )

    emit:
    bed    = BEDTOOLS_SORT.out.sorted       // channel: [ val(meta), path(bed) ]
    bigbed = UCSC_BEDTOBIGBED.out.bigbed    // channel: [ val(meta), path(bigBed) ]
}
