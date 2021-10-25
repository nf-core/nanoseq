// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SAMTOOLS_VIEW_BAM {
    tag "$meta.id"
    label 'process_medium'
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::samtools=1.10" : null)
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(sam)
    
    output:
    tuple val(meta), path("*.bam") ,emit: bam

    script:
    """
    samtools view -b -h -O BAM -@ $task.cpus -o ${meta.id}.bam $sam
    """
}
