/*
 * Alignment with GRAPHMAP2
 */

include { GRAPHMAP2_INDEX    } from '../../../modules/nf-core/graphmap2/index/main'
include { GRAPHMAP2_ALIGN    } from '../../../modules/nf-core/graphmap2/align/main'
include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow ALIGN_GRAPHMAP2 {
    take:
    ch_fasta_fai // channel (value): [ val(meta), path(fasta), path(fai) ]
    ch_fastq     // channel: [ val(meta), path(fastq) ]
    ch_gtf       // channel (value): path(gtf) for spliced alignment, or [] for unspliced

    main:
    ch_fasta = ch_fasta_fai.map { _meta, fasta, _fai -> fasta }

    /*
     * Create genome/transcriptome index
     */
    GRAPHMAP2_INDEX(ch_fasta, ch_gtf)

    /*
     * Map reads with GRAPHMAP2 (SAM output)
     */
    GRAPHMAP2_ALIGN(ch_fastq, ch_fasta, GRAPHMAP2_INDEX.out.index, ch_gtf)

    /*
     * Coordinate-sort SAM into BAM and write the BAI index in the same step
     */
    SAMTOOLS_SORT(GRAPHMAP2_ALIGN.out.sam, ch_fasta_fai, 'bai')

    /*
     * SAMtools stats, flagstat and idxstats
     */
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.index, by: 0)
    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    emit:
    index    = GRAPHMAP2_INDEX.out.index         // channel: path(gmidx)
    bam      = SAMTOOLS_SORT.out.bam             // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_SORT.out.index           // channel: [ val(meta), path(bai) ]
    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]
}
