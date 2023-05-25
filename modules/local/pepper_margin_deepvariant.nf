process PEPPER_MARGIN_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    if (params.deepvariant_gpu) {
        container 'docker.io/kishwars/pepper_deepvariant:r0.8-gpu'
    } else {
        container 'docker.io/kishwars/pepper_deepvariant:r0.8'
    }

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("*vcf.gz")    ,  emit: vcf
    tuple val(meta), path("*vcf.gz.tbi"),  emit: tbi
    path "versions.yml"                 ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix      = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args ?: ""
    def gpu     = params.deepvariant_gpu ? "-g" : ""
    def regions = intervals ? "--regions ${intervals}" : ""

    """
    mkdir -p "${prefix}"

    run_pepper_margin_deepvariant call_variant \\
        -b "${input}" \\
        -f "${fasta}" \\
        -o "." \\
        -p "${prefix}" \\
        -t ${task.cpus} \\
        $gpu \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_margin_deepvariant: \$(run_pepper_margin_deepvariant --version | sed 's/VERSION: //g')
    END_VERSIONS
    """
}
