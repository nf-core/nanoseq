/*
 * RNA MODIFICATION DETECTION WITH XPORE AND M6ANET
 */

include { NANOPOLISH_INDEX_EVENTALIGN } from '../../modules/local/nanopolish_index_eventalign'
include { XPORE_DATAPREP        } from '../../modules/local/xpore_dataprep'
include { XPORE_DIFFMOD         } from '../../modules/local/xpore_diffmod'
include { M6ANET_DATAPREP        } from '../../modules/local/m6anet_dataprep'
include { M6ANET_INFERENCE      } from '../../modules/local/m6anet_inference'

workflow RNA_MODIFICATION_XPORE_M6ANET {
    take:
    ch_nanopolish_bam_fast5

    main:

    /*
     * Align current signals to reference with Nanopolish
     */
    NANOPOLISH_INDEX_EVENTALIGN { ch_nanopolish_bam_fast5 }
    ch_nanopolish_outputs = NANOPOLISH_INDEX_EVENTALIGN.out.nanopolish_outputs
    nanopolish_version    = NANOPOLISH_INDEX_EVENTALIGN.out.versions

    xpore_version = ''
    ch_xpore_dataprep_dirs = ''
    if (!params.skip_xpore) {

        /*
         * Prepare data with xpore
         */
        XPORE_DATAPREP( ch_nanopolish_outputs )
        ch_xpore_dataprep_dirs = XPORE_DATAPREP.out.dataprep_outputs
        ch_xpore_dataprep_dirs
            .map{ it -> it[1]+','+it[0].id }
            .set{ ch_xpore_diffmod_inputs }

        /*
         * Differential modification expression with xpore
         */
        XPORE_DIFFMOD{ ch_xpore_diffmod_inputs.collect() }
        xpore_version    = XPORE_DIFFMOD.out.versions
    }

    /*
    * Detect m6A sites with m6anet
    */
    m6anet_version = ''
    if (!params.skip_m6anet) {
        /*
        * Detect m6A sites with m6anet
        */
        M6ANET_DATAPREP( ch_nanopolish_outputs )
        ch_m6anet_dataprep_dirs = M6ANET_DATAPREP.out.dataprep_outputs
        M6ANET_INFERENCE{ ch_m6anet_dataprep_dirs }
        m6anet_version   = M6ANET_INFERENCE.out.versions
    }

    emit:
    ch_nanopolish_outputs
    nanopolish_version
    ch_xpore_dataprep_dirs
    xpore_version
    m6anet_version
}
