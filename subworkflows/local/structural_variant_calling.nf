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

    ch_versions         = Channel.empty()

    /*
    * Get chromosomes from bam file for splitting calling
    */
    GET_CHROM_NAMES( ch_view_sortbam )
    ch_chrom_names = GET_CHROM_NAMES.out.chrom_names

    /*
     * Map chromosome names
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
     * Call structural variants with sniffles
     */
    if (params.structural_variant_caller == 'sniffles') {

        /*
         * Call structural variants with sniffles
         */
        SNIFFLES( ch_view_sortbam )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions)

        /*
         * Zip sniffles vcf
         */
        SNIFFLES_BGZIP_VCF( SNIFFLES.out.sv_calls )
        ch_sv_calls_vcf = SNIFFLES_BGZIP_VCF.out.gz
        ch_versions = ch_versions.mix(SNIFFLES_BGZIP_VCF.out.versions)

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
        CUTESV( ch_view_sortbam, ch_fasta )
        ch_versions = ch_versions.mix(CUTESV.out.versions)

        /*
         * Zip cutesv vcf
         */
        CUTESV_BGZIP_VCF( CUTESV.out.sv_calls )
        ch_sv_calls_vcf = CUTESV_BGZIP_VCF.out.gz
        ch_versions = ch_versions.mix(CUTESV_BGZIP_VCF.out.versions)

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
