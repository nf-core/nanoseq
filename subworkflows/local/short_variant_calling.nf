/*
 * Short variant calling test
 */

include { CLAIR3                                } from '../../modules/local/clair3'
include { DEEPVARIANT                           } from '../../modules/nf-core/deepvariant/main'
include { PEPPER_MARGIN_DEEPVARIANT             } from '../../modules/local/pepper_margin_deepvariant'

workflow SHORT_VARIANT_CALLING {

    take:
    ch_sorted_bam // channel: [mandatory] [ meta, bam ]
    ch_sorted_bai // channel: [mandatory] [ meta, bai ]
    ch_fasta      // channel: [mandatory] [ meta, fasta ]
    ch_fai        // channel: [mandatory] [ meta, fai ]

    main:
    ch_short_calls_vcf              = Channel.empty()
    ch_short_calls_vcf_tbi          = Channel.empty()
    ch_short_calls_gvcf             = Channel.empty()
    ch_short_calls_gvcf_tbi         = Channel.empty()
    ch_versions                     = Channel.empty()

    // Map inputs
    ch_sorted_bam
        .join(ch_sorted_bai, by: 0)
        .map{ meta, bam, bai -> [ meta, bam, bai, [] ] }
        .set { ch_input }

    /*
     * Call short variants
     */
    if (params.variant_caller == 'clair3') {

        /*
         * Call short variants with clair3
         */
        CLAIR3 ( ch_input, ch_fasta, ch_fai, [[ id:'null' ],[]] )
        ch_short_calls_vcf = CLAIR3.out.vcf
        ch_short_calls_vcf_tbi = CLAIR3.out.tbi
        ch_versions = ch_versions.mix(CLAIR3.out.versions)

    } else if (params.variant_caller == 'deepvariant') {

        /*
         * Call variants with deepvariant
         */
        DEEPVARIANT( ch_input, ch_fasta, ch_fai, [[ id:'null' ],[]] )
        ch_short_calls_vcf  = DEEPVARIANT.out.vcf
        ch_short_calls_vcf  = DEEPVARIANT.out.vcf_tbi
        ch_short_calls_gvcf = DEEPVARIANT.out.gvcf
        ch_short_calls_gvcf = DEEPVARIANT.out.gvcf_tbi
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)


    } else {

        /*
         * Call variants with pepper_margin_deepvariant
         */
        PEPPER_MARGIN_DEEPVARIANT( ch_input, ch_fasta, ch_fai, [[ id:'null' ],[]] )
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
