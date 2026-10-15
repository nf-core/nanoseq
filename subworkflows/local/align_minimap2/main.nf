/*
 * Alignment with MINIMAP2
 */

include { MINIMAP2_INDEX     } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN     } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow ALIGN_MINIMAP2 {
    take:
    ch_fasta_fai // channel (value): [ val(meta), path(fasta), path(fai) ]
    ch_fastq     // channel: [ val(meta), path(fastq) ]

    main:
    /*
     * Create genome/transcriptome index
     */
    MINIMAP2_INDEX(ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] })

    /*
     * Map reads with MINIMAP2; the module pipes into 'samtools sort --write-index',
     * so a coordinate-sorted BAM and its BAI are produced in a single step
     */
    MINIMAP2_ALIGN(
        ch_fastq,
        MINIMAP2_INDEX.out.index,
        true,  // bam_format
        'bai', // bam_index_extension
        false, // cigar_paf_format
        false, // cigar_bam
    )

    /*
     * SAMtools stats, flagstat and idxstats
     */
    ch_bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, by: 0)
    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    emit:
    index    = MINIMAP2_INDEX.out.index          // channel: [ val(meta), path(mmi) ]
    bam      = MINIMAP2_ALIGN.out.bam            // channel: [ val(meta), path(bam) ]
    bai      = MINIMAP2_ALIGN.out.index          // channel: [ val(meta), path(bai) ]
    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]
}
