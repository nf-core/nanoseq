// Import generic module functions
include { saveFiles; getProcessName } from './functions'

process BAM_RENAME {
    tag "$meta.id"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"
    } else {
        container "quay.io/baselibrary/ubuntu:latest"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    """
    [ ! -f ${meta.id}.bam ] && ln -s $bam ${meta.id}.bam
    """
}
