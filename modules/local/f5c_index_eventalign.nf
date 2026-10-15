process F5C_INDEX_EVENTALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::f5c=1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/f5c:1.5--h56e2c18_1' :
        'quay.io/biocontainers/f5c:1.5--h56e2c18_1' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fastqgz), path(bam), path(bai)

    output:
    tuple val(meta), path(genome), path(gtf), path("*eventalign.txt"), path("*summary.txt"), emit: f5c_outputs
    tuple val("${task.process}"), val('f5c'), eval("f5c --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_f5c, topic: versions

    script:
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def sample_summary    = "${prefix}_summary.txt"
    def sample_eventalign = "${prefix}_eventalign.txt"
    def fastq             = "${prefix}.fastq"
    def fast5             = "${meta.fast5}"
    """
    gunzip -c $fastqgz > $fastq
    f5c index -d $fast5 $fastq
    f5c eventalign \\
        --reads $fastq \\
        --bam $bam \\
        --genome $genome \\
        --scale-events \\
        --signal-index \\
        --rna \\
        --min-mapq 0 \\
        --min-recalib-events 1 \\
        --summary $sample_summary \\
        --threads $task.cpus \\
        > $sample_eventalign
    """
}
