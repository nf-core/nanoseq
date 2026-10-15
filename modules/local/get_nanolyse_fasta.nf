process GET_NANOLYSE_FASTA {
    label "process_single"

    conda "conda-forge::curl=7.88.1"
    // Debian bookworm image with bash, git (2.39.5) and curl (7.88.1)
    container "docker.io/library/buildpack-deps:bookworm-scm@sha256:8ea58fbc817960cd62d512044956ff00cb94865e549eeec32c12139413262004"

    output:
    path "*fasta.gz", emit: fasta
    tuple val("${task.process}"), val('curl'), eval("curl --version | sed '1!d;s/^curl //;s/ .*//'"), emit: versions_curl, topic: versions

    script:
    """
    curl \\
        -L https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz \\
        -o lambda.fasta.gz
    """
}
