/*
 * RNA MODIFICATION DETECTION WITH XPORE AND M6ANET
 */

include { NANOPOLISH_EVENTALIGN } from '../../../modules/local/nanopolish/eventalign/main'
include { XPORE_DATAPREP              } from '../../../modules/local/xpore/dataprep/main'
include { XPORE_DIFFMOD               } from '../../../modules/local/xpore/diffmod/main'
include { M6ANET_DATAPREP             } from '../../../modules/local/m6anet/dataprep/main'
include { M6ANET_INFERENCE            } from '../../../modules/local/m6anet/inference/main'

workflow RNA_MODIFICATION_XPORE_M6ANET {
    take:
    ch_nanopolish_input // channel: [ val(meta), path(genome), path(gtf), path(fastq), path(bam), path(bai), path(fast5_dir) ]
    skip_xpore          // val: boolean
    skip_m6anet         // val: boolean

    main:
    /*
     * Align current signals to reference with Nanopolish
     */
    NANOPOLISH_EVENTALIGN( ch_nanopolish_input )
    ch_nanopolish_outputs = NANOPOLISH_EVENTALIGN.out.nanopolish_outputs

    ch_xpore_dataprep = channel.empty()
    ch_xpore_diffmod  = channel.empty()
    if (!skip_xpore) {
        /*
         * Prepare data with xpore
         */
        XPORE_DATAPREP( ch_nanopolish_outputs )
        ch_xpore_dataprep = XPORE_DATAPREP.out.dataprep_outputs

        /*
         * Differential modification expression with xpore
         * Input: [ [ condition, ... ], [ dataprep_dir, ... ] ] (one condition per sample)
         */
        ch_xpore_diffmod_input = ch_xpore_dataprep
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { samples ->
                [ samples.collect { sample -> sample[0].id }, samples.collect { sample -> sample[1] } ]
            }
        XPORE_DIFFMOD( ch_xpore_diffmod_input )
        ch_xpore_diffmod = XPORE_DIFFMOD.out.diffmod_outputs
    }

    ch_m6anet_outputs = channel.empty()
    if (!skip_m6anet) {
        /*
         * Detect m6A sites with m6anet
         */
        M6ANET_DATAPREP( ch_nanopolish_outputs )
        M6ANET_INFERENCE( M6ANET_DATAPREP.out.dataprep_outputs )
        ch_m6anet_outputs = M6ANET_INFERENCE.out.m6anet_outputs
    }

    emit:
    nanopolish_outputs = ch_nanopolish_outputs // channel: [ val(meta), path(genome), path(gtf), path(eventalign), path(summary) ]
    xpore_dataprep     = ch_xpore_dataprep     // channel: [ val(meta), path(dataprep_dir) ]
    xpore_diffmod      = ch_xpore_diffmod      // channel: path(diffmod*)
    m6anet_outputs     = ch_m6anet_outputs     // channel: path(results_dir)
}
