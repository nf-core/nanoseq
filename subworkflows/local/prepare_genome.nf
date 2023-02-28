/*
 * Prepare genome/transcriptome before alignment
 */
include { CUSTOM_GETCHROMSIZES } from '../../modules/nf-core/custom/getchromsizes/main'
include { GTF2BED              } from '../../modules/local/gtf2bed'

workflow PREPARE_GENOME {

    main:

    ch_versions = Channel.empty()
    ch_fasta    = Channel.empty()
    ch_fai      = Channel.empty()
    ch_sizes    = Channel.empty()
    ch_gtf      = Channel.empty()
    ch_bed      = Channel.empty()

    if (params.fasta) {

        Channel.fromPath(params.fasta)
            .collect()
            .set { ch_fasta }

        /*
         * Generates a FASTA file of chromosome sizes and a fasta index file
         */
        CUSTOM_GETCHROMSIZES (ch_fasta)
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.collect()
        ch_fai = CUSTOM_GETCHROMSIZES.out.fai.collect()
        ch_versions = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    }

    if (params.gtf) {

        Channel.fromPath(params.gtf)
            .collect()
            .set { ch_gtf }

        /*
         * Convert GTF to BED12
         */
        GTF2BED (ch_gtf)
        ch_bed = GTF2BED.out.gtf_bed.collect()
        ch_versions = ch_versions.mix(GTF2BED.out.versions)

    }

    emit:
    ch_fasta
    ch_fai
    ch_sizes
    ch_gtf
    ch_bed
    ch_versions
}
