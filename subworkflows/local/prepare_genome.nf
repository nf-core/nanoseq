/*
 * Prepare genome/transcriptome before alignment
 */

include { GET_CHROM_SIZES  } from '../../modules/local/get_chrom_sizes'
include { GTF2BED          } from '../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX   } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    fasta
    gtf

    main:

    /*
     * Make chromosome sizes file
     */
    if (params.fasta) {
        GET_CHROM_SIZES (fasta)
        ch_chrom_sizes = GET_CHROM_SIZES.out.sizes
        samtools_version = GET_CHROM_SIZES.out.versions
        ch_chrom_sizes.view()

    /*
     * Make fasta index
     */
    //SAMTOOLS_FAIDX (fasta)
    //ch_fai = SAMTOOLS_FAIDX.out.fai

    /*
     * Convert GTF to BED12
     */
    //GTF2BED (gtf)
    //ch_gtf_bed = GTF2BED.out.gtf_bed
    //gtf2bed_version = GTF2BED.out.versions

    emit:
    ch_chrom_sizes
    //ch_fai
    //ch_gtf_bed
    //samtools_version
    //gtf2bed_version
}
