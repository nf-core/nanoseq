process MEDAKA_VARIANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::medaka=1.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0' :
        'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)

    output:
    tuple val(meta), path ("$output_vcf"), emit: vcf // vcf files
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def split_mnps       =  params.split_mnps   ? "-l"                        : ''
    def phase_vcf        =  params.phase_vcf    ? "-p"                        : ''

    output_dir = "${meta.id}"
    output_vcf = output_dir+"/round_1.vcf"
    """

    medaka_variant \\
        -d \\
        -f $fasta \\
        -i $input \\
        -o $output_dir \\
        -t $task.cpus \\
        $split_mnps \\
        $phase_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
