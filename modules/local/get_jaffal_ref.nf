process GET_JAFFAL_REF {
    label "process_single"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    output:
    tuple val(null), path("for_jaffal.tar.gz"), emit: ch_jaffal_ref
    path "versions.yml"                       , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    curl \\
    -L https://ndownloader.figshare.com/files/28168755 \\
    -o for_jaffal.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | grep "curl" | sed "s/curl //; s/ .*\$//")
    END_VERSIONS
    """
}
