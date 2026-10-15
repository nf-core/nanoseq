process NANOPOLISH_INDEX_EVENTALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanopolish=0.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.13.2--he3b7ca5_2' :
        'quay.io/biocontainers/nanopolish:0.13.2--he3b7ca5_2' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fastq), path(bam), path(bai), path(fast5, stageAs: 'fast5_dir')

    output:
    tuple val(meta), path(genome), path(gtf), path("*eventalign.txt"), path("*summary.txt"), emit: nanopolish_outputs
    tuple val("${task.process}"), val('nanopolish'), eval("nanopolish --version | sed '1!d;s/.* //'"), emit: versions_nanopolish, topic: versions

    script:
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def sample_summary    = "${prefix}_summary.txt"
    def sample_eventalign = "${prefix}_eventalign.txt"
    """
    nanopolish index -d $fast5 $fastq
    nanopolish eventalign \\
        --reads $fastq \\
        --bam $bam \\
        --genome $genome \\
        --scale-events \\
        --signal-index \\
        --summary $sample_summary \\
        --threads $task.cpus \\
        > $sample_eventalign
    """

    stub:
    def prefix            = task.ext.prefix ?: "${meta.id}"
    def sample_summary    = "${prefix}_summary.txt"
    def sample_eventalign = "${prefix}_eventalign.txt"
    """
    touch $sample_summary
    touch $sample_eventalign
    """
}
