/*
 * Short variant calling test
 */

include { CLAIR3                        } from '../../modules/local/clair3'
include { TABIX_BGZIP as CLAIR3_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as CLAIR3_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'
include { DEEPVARIANT } from '../../modules/nf-core/deepvariant/main'
include { TABIX_TABIX as DEEPVARIANT_TABIX_VCF  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as DEEPVARIANT_TABIX_GVCF } from '../../modules/nf-core/tabix/tabix/main'
include { PEPPER_MARGIN_DEEPVARIANT             } from '../../modules/local/pepper_margin_deepvariant'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_sorted_bam
    ch_sorted_bai
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

    ch_sorted_bam
         .join(ch_sorted_bai, by: 0)
         .map { it -> [ it[0], it[1], it[2], [] ] }
         .view()
         .set { ch_shortv_input }
    ch_sorted_bam
         .combine(ch_fasta.map{it->it[1]})
         .map { it -> it[2] }
         .set { ch_fasta }
    ch_sorted_bam
         .combine(ch_fai.map{it->it[1]})
         .map { it -> it[2] }
         .set { ch_fai }

    if (params.variant_caller == 'clair3') {

        /*
         * Call short variants with medaka
         */
        CLAIR3 ( ch_shortv_input.map{ it -> [ it[0], it[1], it[2] ] }, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix(medaka_version = CLAIR3.out.versions)

        /*
         * Zip medaka vcf
         */
        CLAIR3_BGZIP_VCF( CLAIR3.out.vcf )
        ch_short_calls_vcf  = CLAIR3_BGZIP_VCF.out.output
        ch_versions = ch_versions.mix(bgzip_version = CLAIR3_BGZIP_VCF.out.versions)

        /*
         * Index medaka vcf.gz
         */
        CLAIR3_TABIX_VCF( ch_short_calls_vcf )
        ch_short_calls_vcf_tbi  = CLAIR3_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(tabix_version = CLAIR3_TABIX_VCF.out.versions)

    } else if (params.variant_caller == 'deepvariant') {

        /*
         * Call variants with deepvariant
         */
        ch_sorted_bam
             .join(ch_sorted_bai, by: 0)
             .map { it -> [ it[0], it[1], it[2], [] ] }
             .view()
             .set { ch_deepvariant_input }
        ch_sorted_bam
             .combine(ch_fasta.map{it->it[1]})
             .map { it -> it[2] }
             .set { ch_fasta }
        ch_sorted_bam
             .combine(ch_fai.map{it->it[1]})
             .map { it -> it[2] }
             .set { ch_fai }

        DEEPVARIANT( ch_deepvariant_input, ch_fasta, ch_fai )
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
        PEPPER_MARGIN_DEEPVARIANT( ch_shortv_input, ch_fasta, ch_fai )
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
