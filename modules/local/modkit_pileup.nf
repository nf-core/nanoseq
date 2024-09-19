process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_high_memory'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-modkit:0.4.0--h5c23e0d_0' :
        'quay.io/biocontainers/ont-modkit:0.4.0--h5c23e0d_0' }"

    input:
    tuple val(meta), path(aligned_mod_bam), path(bai)

    output:
    tuple val(meta), path(bedmethyl)     , emit: bedmethyl
    path "versions.yml"                  , emit: versions

    script:
    bedmethyl = "$meta.id" +".bed"
    """
    modkit pileup $aligned_mod_bam $bedmethyl --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(echo \$(modkit --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS

    gzip basecall.fastq
    """
}

