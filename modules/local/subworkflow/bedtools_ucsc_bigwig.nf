/*
 * Convert BAM to BigWig
 */

params.bigwig_options   = [:]

include { BEDTOOLS_GENOMECOV    } from '../process/bedtools_genomecov'                     addParams( options: params.bigwig_options )
include { UCSC_BEDGRAPHTOBIGWIG } from '../../nf-core/software/ucsc/bedgraphtobigwig/main' addParams( options: params.bigwig_options )

workflow BEDTOOLS_UCSC_BIGWIG {
    take:
    ch_sortbam // channel: [ val(meta), [ reads ] ]
    
    main:
    /*
     * Convert BAM to BEDGraph
     */
     BEDTOOLS_GENOMECOV ( ch_sortbam )
     ch_bedgraph = BEDTOOLS_GENOMECOV.out.bedgraph

    /*
     * Convert BEDGraph to BigWig
     */

    UCSC_BEDGRAPHTOBIGWIG ( ch_bedgraph )

    emit:
    bigwig      = UCSC_BEDGRAPHTOBIGWIG.bigWig
}
