process GET_JAFFAL_REF {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    output:
    path "for_jaffal.tar.gz"  , emit: ch_jaffal_ref

    script:
    """
    curl \\
    -L https://ndownloader.figshare.com/files/28168755 \\
    -o for_jaffal.tar.gz
    """
}
