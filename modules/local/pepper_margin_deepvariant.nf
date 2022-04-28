process PEPPER_MARGIN_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/kishwars/pepper_deepvariant:r0.8' :
        'docker.io/kishwars/pepper_deepvariant:r0.8' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${prefix}/*vcf.gz")     ,  emit: vcf
    tuple val(meta), path("${prefix}/*vcf.gz.tbi") ,  emit: tbi
    path "versions.yml"                            ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ""
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""
    //def gvcf    = params.make_gvcf ? "--gvcf" : ""

    """
    mkdir -p "${prefix}"

    run_pepper_margin_deepvariant call_variant \\
        -b "${input}" \\
        -f "${fasta}" \\
        -o "${prefix}" \\
        -p "${prefix}" \\
        -t ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_margin_deepvariant: \$(run_pepper_margin_deepvariant --version | sed 's/VERSION: //g')
    END_VERSIONS
    """
}
