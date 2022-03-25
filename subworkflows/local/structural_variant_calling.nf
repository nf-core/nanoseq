/*
 * Structural variant calling
 */

include { GET_CHROM_NAMES                       } from '../../modules/local/get_chrom_names'
include { SNIFFLES                              } from '../../modules/local/sniffles'
include { TABIX_BGZIP as SNIFFLES_BGZIP_VCF     } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_TABIX as SNIFFLES_TABIX_VCF     } from '../../modules/nf-core/modules/tabix/tabix/main'
include { CUTESV                                } from '../../modules/local/cutesv'
include { TABIX_BGZIP as CUTESV_BGZIP_VCF       } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_TABIX as CUTESV_TABIX_VCF       } from '../../modules/nf-core/modules/tabix/tabix/main'


workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai

    main:
    ch_sv_calls_vcf     = Channel.empty()
    ch_sv_calls_vcf_tbi = Channel.empty()

    sniffles_version    = Channel.empty()
    cutesv_version      = Channel.empty()
    bgzip_version       = Channel.empty()
    tabix_version       = Channel.empty()

    /*
    * Get names of chromosomes from bam file
    */
    GET_CHROM_NAMES( ch_view_sortbam )
    ch_chrom_names = GET_CHROM_NAMES.out.chrom_names

    /*
    * Map
    */
    ch_chrom_names
    .splitCsv()
    .combine( ch_view_sortbam )
    .unique()
    .map { chroms, meta, sizes, is_transcripts, bam, bai ->
    new_meta = meta.clone()
    new_meta.id = meta.id + "_" + chroms
    new_meta.sample = meta.id
    new_meta.chr = chroms
    [new_meta, bam, bai, chroms]
    }.set { ch_view_sortbam_split }

    /*
    * Call structural variants
    */
    if (params.structural_variant_caller == 'sniffles') {
        /*
        * SNIFFLES
        */
        SNIFFLES( ch_view_sortbam )
        sniffles_version = SNIFFLES.out.versions

        SNIFFLES_BGZIP_VCF( SNIFFLES.out.sv_calls )
        ch_sv_calls_vcf = SNIFFLES_BGZIP_VCF.out.gz
        bgzip_version = SNIFFLES_BGZIP_VCF.out.versions

        SNIFFLES_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = SNIFFLES_TABIX_VCF.out.tbi
        tabix_version = SNIFFLES_TABIX_VCF.out.versions

    } else {
        /*
        * CUTESV
        */
        CUTESV( ch_view_sortbam, ch_fasta )
        cutesv_version = CUTESV.out.versions

        CUTESV_BGZIP_VCF( CUTESV.out.sv_calls )
        ch_sv_calls_vcf = CUTESV_BGZIP_VCF.out.gz
        bgzip_version = CUTESV_BGZIP_VCF.out.versions

        CUTESV_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = CUTESV_TABIX_VCF.out.tbi
        tabix_version = CUTESV_TABIX_VCF.out.versions

    }

    emit:
    ch_sv_calls_vcf
    ch_sv_calls_vcf_tbi
    bgzip_version
    tabix_version
    sniffles_version
    cutesv_version

}
