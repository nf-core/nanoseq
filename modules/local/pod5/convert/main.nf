process FAST5_TO_POD5 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a15a8b329ccf4b1ecaeb3695139ce4b9819573cb2b8c836d858fd76cb51c8f2e/data' :
        'community.wave.seqera.io/library/pod5:0.3.47--e5110e3814a7ce5c' }"

    input:
    tuple val(meta), path(input_path)

    output:
    tuple val(meta), path("pod5/*"), emit: pod5
    tuple val("${task.process}"), val('pod5'), eval("pod5 --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_pod5, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pod5 convert fast5 \\
        $args \\
        --threads $task.cpus \\
        $input_path \\
        --output pod5/${prefix}.pod5
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p pod5
    touch pod5/${prefix}.pod5
    """
}
