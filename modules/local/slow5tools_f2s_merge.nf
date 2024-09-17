process SLOW5_F2S_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::slow5tools==1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/slow5tools:1.2.0--h56e2c18_1' :
        'quay.io/biocontainers/slow5tools:1.2.0--h56e2c18_1' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fastq), path(bam), path(bai)

    output:
    tuple val(meta), path(genome), path(gtf), path("*eventalign.txt"), path("*summary.txt"), emit: nanopolish_outputs
    path "versions.yml"                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    slow5tools f2s fast5_dir -d blow5_dir
    slow5tools merge blow5_dir -o file.blow5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$( slow5tools -V | head -n 1 | tail -c 6 )
    END_VERSIONS
    """
}
