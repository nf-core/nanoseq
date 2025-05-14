/*
 * RNA FUSION DETECTION WITH JAFFAL
 */

include { GET_JAFFAL_REF } from '../../modules/local/get_jaffal_ref'
include { UNTAR          } from '../../modules/nf-core/untar/main'
include { JAFFAL         } from '../../modules/local/jaffal'

workflow RNA_FUSIONS_JAFFAL {
    take:
    ch_sample_simple
    jaffal_ref_dir

    main:
    if (jaffal_ref_dir) {
        ch_jaffal_ref_dir = file(params.jaffal_ref_dir, checkIfExists: true)
    } else {

        /*
         * Get jaffel reference
         */
        GET_JAFFAL_REF()

        /*
         * Untar jaffel reference
         */
        UNTAR( GET_JAFFAL_REF.out.ch_jaffal_ref )
        UNTAR.out.untar
            .map { it -> [ it[1] ]}
            .set { ch_jaffal_ref_dir }

    ch_jaffal_ref_dir.view()

    }

    ch_sample_simple.view()

    /*
    * Align current signals to reference with jaffel
    */
    JAFFAL( ch_sample_simple, ch_jaffal_ref_dir )
}
