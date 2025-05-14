process BLUE_CRAB {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::slow5tools==1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/slow5tools:1.2.0--h56e2c18_1' :
        'quay.io/biocontainers/slow5tools:1.2.0--h56e2c18_1' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fastq), path(bam), path(bai), path(pod5)

    output:
    tuple val(meta), path(genome), path(gtf), path(fastq), path(bam), path(bai), path(blow5), emit: nanopolish_outputs
    path "versions.yml"                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    blue-crab p2s $pod5 -o $blow5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blue-crab: \$( blue-crab -V | tail -c 6 )
    END_VERSIONS
    """
}
