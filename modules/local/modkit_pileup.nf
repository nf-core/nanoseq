process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high_memory'

    conda "bioconda::ont-modkit=0.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.4.0--h5c23e0d_0' :
        'quay.io/biocontainers/ont-modkit:0.4.0--h5c23e0d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bedmethyl
    tuple val("${task.process}"), val('modkit'), eval("modkit --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_modkit, topic: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    modkit pileup \\
        $args \\
        $bam \\
        "${prefix}.bed" \\
        --threads $task.cpus
    """
}
