process DORADO_ALIGNER {
    tag "$meta.id"
    label 'process_medium'

    // No conda package available for dorado
    container "docker.io/ontresearch/dorado:sha9809639e07a927bcc0f584dadd5e59674cf59f3f"

    input:
    tuple val(meta), path(reads)
    path fasta

    output:
    tuple val(meta), path("${prefix}.bam")    , emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai
    tuple val("${task.process}"), val('dorado'), eval("dorado --version 2>/dev/null | sed '1!d;s/[+].*//'"), emit: versions_dorado, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    script:
    def args = task.ext.args ?: '--mm2-preset map-ont'
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    dorado aligner \\
        $args \\
        $fasta \\
        $reads \\
        | samtools sort -@ $task.cpus -o ${prefix}.bam -

    samtools index ${prefix}.bam
    """
}
