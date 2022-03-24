/*
 * Short variant calling test
 */

include { GET_CHROM_NAMES                       } from '../../modules/local/get_chrom_names'
include { MEDAKA_VARIANT                        } from '../../modules/local/medaka_variant'
include { TABIX_BGZIP as MEDAKA_BGZIP_VCF       } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { TABIX_TABIX as MEDAKA_TABIX_VCF       } from '../../modules/nf-core/modules/tabix/tabix/main'
include { DEEPVARIANT                           } from '../../modules/local/deepvariant'
include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF  } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF } from '../../modules/nf-core/modules/tabix/tabix/main'
include { SNIFFLES                              } from '../../modules/local/sniffles'
include { CUTESV                                } from '../../modules/local/cutesv'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai

    main:
    ch_short_calls_vcf      = Channel.empty()
    ch_short_calls_vcf_tbi  = Channel.empty()
    ch_short_calls_gvcf     = Channel.empty()
    ch_short_calls_gvcf_tbi = Channel.empty()

    medaka_version          = Channel.empty()
    deepvariant_version     = Channel.empty()
    bgzip_version           = Channel.empty()
    tabix_version           = Channel.empty()

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
    
    //ch_view_sortbam
    //.combine( ch_chrom_names )
    //.unique()
    //.map{ meta, sizes, is_transcripts, bam, bai, chroms ->
    //[meta, bam, bai, chroms ]
    //}.set{ ch_view_sortbam_split }

    /*
    * Call short variants
    */
    if (params.variant_caller == 'medaka'){
        /*
        * MEDAKA
        */
        MEDAKA_VARIANT( ch_view_sortbam_split, ch_fasta )
        medaka_version = MEDAKA_VARIANT.out.versions

        MEDAKA_BGZIP_VCF( MEDAKA_VARIANT.out.vcf )
        ch_short_calls_vcf  = MEDAKA_BGZIP_VCF.out.gz
        bgzip_version = MEDAKA_BGZIP_VCF.out.versions

        MEDAKA_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = MEDAKA_TABIX_VCF.out.tbi
        tabix_version = MEDAKA_TABIX_VCF.out.versions

    } else {
        /*
        * DEEPVARIANT
        */
        DEEPVARIANT( ch_view_sortbam_split, ch_fasta, ch_fai )
        ch_short_calls_vcf  = DEEPVARIANT.out.vcf
        ch_short_calls_gvcf = DEEPVARIANT.out.gvcf
        deepvariant_version = DEEPVARIANT.out.versions

        DEEPVARIANT_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = DEEPVARIANT_TABIX_VCF.out.tbi
        tabix_version = DEEPVARIANT_TABIX_VCF.out.versions

        DEEPVARIANT_TABIX_GVCF( ch_short_calls_gvcf )
        ch_short_calls_gvcf_tbi  = DEEPVARIANT_TABIX_GVCF.out.tbi
        tabix_version = DEEPVARIANT_TABIX_VCF.out.versions

    }


    emit:
    ch_short_calls_vcf
    ch_short_calls_vcf_tbi
    ch_short_calls_gvcf
    ch_short_calls_gvcf_tbi
    bgzip_version
    tabix_version
    medaka_version
    deepvariant_version

}
