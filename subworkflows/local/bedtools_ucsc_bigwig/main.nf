/*
 * Convert BAM to BigWig
 */

include { BEDTOOLS_GENOMECOV    } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDTOOLS_UCSC_BIGWIG {
    take:
    ch_bam   // channel: [ val(meta), path(bam) ]
    ch_sizes // channel (value): [ val(meta), path(sizes) ]

    main:
    ch_chrom_sizes = ch_sizes.map { _meta, sizes -> sizes }

    /*
     * Convert BAM to bedGraph (scale = 1, i.e. no scaling; output sorted with 'sort -k1,1 -k2,2n')
     */
    BEDTOOLS_GENOMECOV(
        ch_bam.map { meta, bam -> [meta, bam, 1] },
        ch_chrom_sizes,
        'bedGraph', // extension
        true,       // sort
    )

    /*
     * Convert bedGraph to BigWig
     */
    UCSC_BEDGRAPHTOBIGWIG(BEDTOOLS_GENOMECOV.out.genomecov, ch_chrom_sizes)

    emit:
    bedgraph = BEDTOOLS_GENOMECOV.out.genomecov  // channel: [ val(meta), path(bedGraph) ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig  // channel: [ val(meta), path(bigWig) ]
}
