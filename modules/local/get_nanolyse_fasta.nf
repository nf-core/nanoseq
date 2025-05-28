process GET_NANOLYSE_FASTA {
    label "process_single"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    output:
    path "*fasta.gz"   , emit: ch_nanolyse_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    curl \\
    -L https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz \\
    -o lambda.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(echo \$(curl --version  2>&1) | grep "curl" | sed "s/curl //; s/ .*\$//")
    END_VERSIONS
    """
}
