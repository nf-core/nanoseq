process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1' :
        'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1' }"

    input:
    tuple val(meta), path(sizes), path(bedgraph)

    output:
    tuple val(meta), path(sizes), path("*.bigWig"), emit: bigwig
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '377'
    """
    bedGraphToBigWig $bedgraph $sizes ${meta.id}.bigWig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc_bedgraphtobigwig: \$(echo $VERSION)
    END_VERSIONS
    """
}
