/*
 * Structural variant calling
 */

include { SNIFFLES                  } from '../../modules/local/sniffles'
include { CUTESV                    } from '../../modules/local/cutesv'

workflow STRUCTURAL_VARIANT_CALLING {

    main:
    take:
    ch_view_sortbam
    ch_index

    main:
    ch_sv_calls = Channel.empty()
    sniffles_version = Channel.empty()
    cutesv_version = Channel.empty()
    if (params.structural_variant_caller == 'sniffles'){
        /*
        * Call structural variants with SNIFFLES
        */
        SNIFFLES( ch_view_sortbam )
        ch_sv_calls = SNIFFLES.out.sv_calls
        sniffles_version = SNIFFLES.out.versions

    } else {
        /*
        * Call structural variants with cuteSV
        */
        CUTESV( ch_view_sortbam, ch_index )
        ch_sv_calls = CUTESV.out.sv_calls
        cutesv_version = CUTESV.out.versions

    }

    emit:
    ch_sv_calls
    sniffles_version
    cutesv_version

}
