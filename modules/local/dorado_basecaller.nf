process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    // No conda package available for dorado
    container "docker.io/ontresearch/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(pod5_path)
    val dorado_device
    val dorado_model
    val dorado_modification

    output:
    tuple val(meta), path("*.fastq.gz"), optional: true, emit: fastq
    tuple val(meta), path("*.bam")     , optional: true, emit: bam
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>/dev/null | sed '1!d;s/[+].*//'"), emit: versions_dorado, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def emit_args = dorado_modification
        ? "--modified-bases ${dorado_modification} > ${prefix}.bam"
        : "--emit-fastq | gzip > ${prefix}.fastq.gz"
    """
    dorado download --model $dorado_model

    dorado basecaller \\
        $args \\
        $dorado_model \\
        $pod5_path \\
        --device $dorado_device \\
        $emit_args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stub_cmd = dorado_modification
        ? "touch ${prefix}.bam"
        : "echo '' | gzip > ${prefix}.fastq.gz"
    """
    $stub_cmd
    """
}
