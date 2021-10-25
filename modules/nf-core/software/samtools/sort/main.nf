// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::samtools=1.10" : null)
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    output_bam = "$meta.id"+".sorted.bam"
    """
    samtools sort $options.args -@ $task.cpus -o $output_bam -T $meta.id $bam
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > samtools.version.txt
    """
}
