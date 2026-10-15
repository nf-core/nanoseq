process GET_TEST_DATA {
    label "process_single"

    conda "conda-forge::git=2.39.2"
    // Debian bookworm image with bash, git (2.39.5) and curl (7.88.1)
    container "docker.io/library/buildpack-deps:bookworm-scm@sha256:8ea58fbc817960cd62d512044956ff00cb94865e549eeec32c12139413262004"

    output:
    path "test-datasets/modification_fast5_fastq/", emit: modification_dir
    tuple val("${task.process}"), val('git'), eval("git --version | sed '1!d;s/.* //'"), emit: versions_git, topic: versions

    script:
    """
    git clone https://github.com/nf-core/test-datasets.git --branch nanoseq --single-branch
    """
}
