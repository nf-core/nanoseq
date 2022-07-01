/*
 * Short variant calling test
 */

include { CLAIR3                                } from '../../modules/local/clair3'
include { DEEPVARIANT                           } from '../../modules/local/deepvariant'
include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF  } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF } from '../../modules/nf-core/modules/tabix/tabix/main'
include { PEPPER_MARGIN_DEEPVARIANT             } from '../../modules/local/pepper_margin_deepvariant'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_view_sortbam
    ch_fasta
    ch_fai

    main:
    ch_short_calls_vcf              = Channel.empty()
    ch_short_calls_vcf_tbi          = Channel.empty()
    ch_short_calls_gvcf             = Channel.empty()
    ch_short_calls_gvcf_tbi         = Channel.empty()
    ch_versions                     = Channel.empty()

    /*
     * Call short variants
     */
    if (params.variant_caller == 'clair3') {

        /*
         * Call short variants with clair3
         */
        CLAIR3( ch_view_sortbam, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix(clair3_version = CLAIR3.out.versions)

    } else if (params.variant_caller == 'deepvariant') {

        /*
        * Call variants with deepvariant
        */
        DEEPVARIANT( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf  = DEEPVARIANT.out.vcf
        ch_short_calls_gvcf = DEEPVARIANT.out.gvcf
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

        /*
         * Index deepvariant vcf.gz
         */
        DEEPVARIANT_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = DEEPVARIANT_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(DEEPVARIANT_TABIX_VCF.out.versions)

        /*
         * Index deepvariant g.vcf.gz
         */
        DEEPVARIANT_TABIX_GVCF( ch_short_calls_gvcf )
        ch_short_calls_gvcf_tbi  = DEEPVARIANT_TABIX_GVCF.out.tbi
        ch_versions = ch_versions.mix(DEEPVARIANT_TABIX_VCF.out.versions)

    } else {

        /*
         * Call variants with pepper_margin_deepvariant (automatic zip + index, docker + singularity only)
         */
        PEPPER_MARGIN_DEEPVARIANT( ch_view_sortbam, ch_fasta, ch_fai )
        ch_short_calls_vcf = PEPPER_MARGIN_DEEPVARIANT.out.vcf
        ch_short_calls_vcf_tbi = PEPPER_MARGIN_DEEPVARIANT.out.tbi
        ch_versions = ch_versions.mix(PEPPER_MARGIN_DEEPVARIANT.out.versions)
    }

    emit:
    ch_short_calls_vcf
    ch_short_calls_vcf_tbi
    ch_short_calls_gvcf
    ch_short_calls_gvcf_tbi
    ch_versions
}
