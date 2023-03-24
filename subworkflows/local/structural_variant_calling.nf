/*
 * Structural variant calling
 */

include { SNIFFLES                              } from '../../modules/nf-core/sniffles/main'
include { BCFTOOLS_SORT as SNIFFLES_SORT_VCF    } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as SNIFFLES_BGZIP_VCF     } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as SNIFFLES_TABIX_VCF     } from '../../modules/nf-core/tabix/tabix/main'
include { CUTESV                                } from '../../modules/nf-core/cutesv/main'
include { BCFTOOLS_SORT as CUTESV_SORT_VCF      } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as CUTESV_BGZIP_VCF       } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as CUTESV_TABIX_VCF       } from '../../modules/nf-core/tabix/tabix/main'


workflow STRUCTURAL_VARIANT_CALLING {

    take:
    ch_sorted_bam
    ch_sorted_bai
    ch_fasta
    ch_fai

    main:
    ch_sv_calls_vcf     = Channel.empty()
    ch_sv_calls_vcf_tbi = Channel.empty()

    ch_versions         = Channel.empty()

    ch_sorted_bam
         .join(ch_sorted_bai, by: 0)
         .map { it -> [ it[0], it[1], it[2] ] }
         .set { ch_sv_input }
    ch_sorted_bam
         .combine(ch_fasta.map{it->it[1]})
         .map { it -> [ [id:'fasta'], it[2] ] }
         .set { ch_fasta }

    /*
     * Call structural variants with sniffles
     */
    if (params.structural_variant_caller == 'sniffles') {

        /*
         * Call structural variants with sniffles
         */
        SNIFFLES( ch_sv_input, ch_fasta )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        SNIFFLES_SORT_VCF( SNIFFLES.out.vcf )
        ch_sv_calls_vcf = SNIFFLES_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(SNIFFLES_SORT_VCF.out.versions)

        /*
         * Index sniffles vcf.gz
         */
        SNIFFLES_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = SNIFFLES_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(SNIFFLES_TABIX_VCF.out.versions)

    } else if (params.structural_variant_caller == 'cutesv') {

        /*
        * Call structural variants with cutesv
        */
        CUTESV( ch_sv_input, ch_fasta )
        ch_versions = ch_versions.mix(CUTESV.out.versions)

        /*
         * Sort structural variants with bcftools
         */
        CUTESV_SORT_VCF( CUTESV.out.vcf )
        ch_sv_calls_vcf = CUTESV_SORT_VCF.out.vcf
        ch_versions = ch_versions.mix(CUTESV_SORT_VCF.out.versions)

        /*
         * Zip cutesv vcf.gz
         */
        CUTESV_TABIX_VCF( ch_sv_calls_vcf )
        ch_sv_calls_tbi  = CUTESV_TABIX_VCF.out.tbi
        ch_versions = ch_versions.mix(CUTESV_TABIX_VCF.out.versions)
    }

    emit:
    ch_sv_calls_vcf
    ch_sv_calls_vcf_tbi
    ch_versions
}

