// Import generic module functions
include { saveFiles; getProcessName } from './functions'

process BAM_RENAME {
    tag "$meta.id"
    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
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
