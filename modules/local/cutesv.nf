process CUTESV {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cutesv=1.0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutesv:1.0.12--pyhdfd78af_0' :
        'quay.io/biocontainers/cutesv:1.0.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)

    output:
    tuple val(meta), path("*_cuteSV.vcf"), emit: sv_calls // vcf files
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cuteSV \
        ${input} \
        ${fasta} \
        ${meta.id}_cuteSV.vcf \
        . \
        --threads $task.cpus \
        --sample ${meta.id} \
        --genotype

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cuteSV: \$( cuteSV --version 2>&1 | sed 's/cuteSV //g' )
    END_VERSIONS
    """
}

