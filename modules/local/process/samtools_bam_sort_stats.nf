// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SAMTOOLS_BAM_SORT_STATS {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::samtools=1.10" : null)
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    input:
    tuple val(sample), path(sizes), val(is_transcripts), path(sam)
    
    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.sorted.bam"), path("*.sorted.bam.bai") ,emit: sortbam
    tuple path("*.flagstat"), path("*.idxstats"), path("*.stats")                                       ,emit: sortbam_stats_multiqc

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -o ${sample}.bam $sam
    samtools sort -@ $task.cpus -o ${sample}.sorted.bam -T $sample ${sample}.bam
    samtools index ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}.sorted.bam.flagstat
    samtools idxstats ${sample}.sorted.bam > ${sample}.sorted.bam.idxstats
    samtools stats ${sample}.sorted.bam > ${sample}.sorted.bam.stats
    """
}
