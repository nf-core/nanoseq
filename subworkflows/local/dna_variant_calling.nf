/*
 * DNA variant calling
 */
params.medaka_vc_options    = [:]
params.sniffles_sv_options    = [:]

include { MEDAKA                  } from '../../modules/local/medaka_variant'      addParams( options: params.medaka_vc_options    )
include { SNIFFLES                  } from '../../modules/local/sniffles'      addParams( options: params.sniffles_sv_options    )

workflow DNA_VARIANT_CALLING {

    main:
    take:
    ch_view_sortbam
    ch_index
    skip_medaka
    skip_sniffles

    main:
    ch_variant_calls = Channel.empty()
    if (!skip_medaka){
        /*
        * Split into a different channel for each chromosome
        */
        // TODO Add module that cuts chromosomes from reference to use as regions to split variant calling 
        //SPLIT_CHROM( ch_view_sortbam )
        //.splitCsv()
        //.combine( ch_view_sortbam ) //
        //.unique()
        //.map { it -> [ it[1], it[4], it[5], it[0] ] } //
        //.map{ meta, bam, bai, chrom -> 
        //new_meta = meta.clone()
        //new_meta.chrom = chrom
        //[new_meta, bam, bai]
        //}.set{ch_bam_vc_chrom}

        /*
        * Call variants with MEDAKA
        */ 
        MEDAKA( ch_view_sortbam, ch_index )
        ch_variant_calls = MEDAKA.out.variant_calls
    }
    ch_sv_calls = Channel.empty()
    if (!skip_sniffles){
        /*
        * Call variants with SNIFFLES
        */
        SNIFFLES( ch_view_sortbam )
        ch_sv_calls = SNIFFLES.out.sv_calls
    } else{
    }

    emit:
    ch_variant_calls
    ch_sv_calls
}