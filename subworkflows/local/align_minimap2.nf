/*
 * Alignment with MINIMAP2
 */

include { MINIMAP2_INDEX as MINIMAP2_INDEX_VARIANT } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_VARIANT } from '../../modules/nf-core/minimap2/align/main'
//include { MINIMAP2_INDEX as MINIMAP2_INDEX_SPLICE  } from '../../modules/nf-core/minimap2/index/main'
//include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SPLICE  } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX                           } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../../subworkflows/nf-core/bam_stats_samtools/main'


workflow ALIGN_MINIMAP2 {
    take:
    ch_fasta // channel: [ path fasta ]
    ch_sizes // channel: [ path sizes ]
    ch_gtf   // channel: [ path gtf ]
    ch_bed   // channel: [ path bed ]
    ch_fastq // channel: [ val(meta), path(fastq) ]

    main:

    ch_versions = Channel.empty()

    if (params.protocol == "DNA") {

        /*
         * Create genome/transcriptome index for variant calling with MINIMAP2
         */
        MINIMAP2_INDEX_VARIANT (ch_fasta)
        ch_index    = MINIMAP2_INDEX_VARIANT.out.index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX_VARIANT.out.versions)

        /*
         * Map reads with MINIMAP2
         */
        MINIMAP2_ALIGN_VARIANT (ch_fastq, ch_fasta, true, false, false)
        ch_bam = MINIMAP2_ALIGN_VARIANT.out.bam // channel: [ val (meta) , path (bam) ] ]
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN_VARIANT.out.versions)
    }

    //if (params.protocol == 'cDNA' || params.protocol == 'directRNA') {

        /*
         * Create genome/transcriptome index for splice variant calling
         */
        //MINIMAP2_INDEX_SPLICE (ch_fasta)
        //ch_index    = MINIMAP2_INDEX_SPLICE.out.index
        //ch_versions = ch_versions.mix(MINIMAP2_INDEX_SPLICE.out.versions)

        /*
         * Map reads with MINIMAP2
         */
        //MINIMAP2_ALIGN (ch_fastq, ch_fasta, ch_index, ch_bed)
        //ch_bam = MINIMAP2_ALIGN.out.bam
        //ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //}

    /*
    * Index mapped reads with SAMTOOLS
    */
    SAMTOOLS_INDEX(ch_bam)
    ch_bai = SAMTOOLS_INDEX.out.bai // channel: [ val (meta) , path (bai) ] ]
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam
        .join(ch_bai, by: 0)
        .set{ ch_bam_bai }

    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta)
    ch_stats    = BAM_STATS_SAMTOOLS.out.stats     // channel: [ val(meta), [ stats ] ]
    ch_flagstat = BAM_STATS_SAMTOOLS.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    ch_idxstats  = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    ch_index
    ch_bam_bai
    ch_stats
    ch_flagstat
    ch_idxstats
    ch_versions
}
