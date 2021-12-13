/*
 * Structural variant calling with SNIFFLES
 */

params.sniffles_sv_options    = [:]

include { SNIFFLES                  } from '../../modules/local/sniffles'      addParams( options: params.sniffles_sv_options    )

workflow SV_SNIFFLES {
    take:
    ch_view_sortbam 

    main:
    /*
     * Call variants with MEDAKA
     */
    SNIFFLES( ch_view_sortbam )
    ch_sv_calls = SNIFFLES.out.sv_calls

    emit:
    ch_sv_calls

}