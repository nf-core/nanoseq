process FAST5_TO_POD5 {
    tag "$meta.id"
    label 'process_medium'

    // No conda directive: pod5 0.2.4 is not available on bioconda/conda-forge
    container "docker.io/yuukiiwa/pod5:0.2.4"

    input:
    tuple val(meta), path(input_path)

    output:
    tuple val(meta), path("pod5/*"), emit: pod5
    tuple val("${task.process}"), val('pod5'), eval("pod5 --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_pod5, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pod5 convert fast5 \\
        $input_path \\
        --output pod5/${prefix}.pod5
    """
}
