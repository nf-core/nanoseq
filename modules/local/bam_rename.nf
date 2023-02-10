process BAM_RENAME {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sed:4.7.0' :
        'quay.io/biocontainers/sed:4.7.0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    [ ! -f ${meta.id}.bam ] && ln -s $bam ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    sed: \$(sed --version 2>&1 | grep "sed (GNU sed)" | sed 's/^.*) //')
    END_VERSIONS
    """
}
