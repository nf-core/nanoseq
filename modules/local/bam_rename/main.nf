process BAM_RENAME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('sed'), eval("sed --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_sed, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -f ${prefix}.bam ]; then
        ln -s $bam ${prefix}.bam
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
