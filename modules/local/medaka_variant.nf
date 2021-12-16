process MEDAKA_VARIANT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::medaka=1.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0' :
        'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(bam), path(bai) //
    path(fasta)

    output:
    path ("${meta.id}/round_1.vcf")    , emit: variant_calls // vcf files
    path "versions.yml"        , emit: versions

    script:
    def args       =  options.args      ?: ''
    def split_mnps =  params.split_mnps ? "-l"                        : ''
    def phase_vcf  =  params.phase_vcf  ? "-p"                        : ''
    """
    medaka_variant \\
        -d \\
        -f $fasta \\
        -i $bam \\
        -o ${meta.id}/ \\
        -t $task.cpus \\
        $split_mnps \\
        $phase_vcf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
