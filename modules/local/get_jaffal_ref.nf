process GET_JAFFAL_REF {
    label "process_single"

    conda "conda-forge::curl=7.88.1"
    // Debian bookworm image with bash, git (2.39.5) and curl (7.88.1)
    container "docker.io/library/buildpack-deps:bookworm-scm@sha256:8ea58fbc817960cd62d512044956ff00cb94865e549eeec32c12139413262004"

    output:
    tuple val([id: 'jaffal_ref']), path("for_jaffal.tar.gz"), emit: ref
    tuple val("${task.process}"), val('curl'), eval("curl --version | sed '1!d;s/^curl //;s/ .*//'"), emit: versions_curl, topic: versions

    script:
    """
    curl \\
        -L https://ndownloader.figshare.com/files/28168755 \\
        -o for_jaffal.tar.gz
    """
}
