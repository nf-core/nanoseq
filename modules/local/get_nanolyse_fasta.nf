process GET_NANOLYSE_FASTA {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    output:
    path "*fasta.gz"  , emit: ch_nanolyse_fasta
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    curl \\
    -L https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz \\
    -o lambda.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | grep "curl" | sed "s/curl //; s/ .*\$//")
    END_VERSIONS
    """
}
