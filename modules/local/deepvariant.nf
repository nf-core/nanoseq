process DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'


    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the DeepVariant tool. Please use docker or singularity containers."
    }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/google/deepvariant:1.3.0' :
        'docker.io/google/deepvariant:1.3.0' }"

    input:
    tuple val(meta), path(input), path(index), val(intervals)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")   ,  emit: vcf
    tuple val(meta), path("${prefix}.g.vcf.gz") ,  emit: gvcf
    path "versions.yml"                         ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions ${intervals}" : ""

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        ${args} \\
        ${regions} \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}
