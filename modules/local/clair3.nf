process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::clair3=0.1.10 conda-forge::python=3.9.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:0.1.10--hdfd78af_0' :
        'biocontainers/clair3:0.1.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${meta.id}/*.gz")    , emit: vcf
    tuple val(meta), path("${meta.id}/*.gz.tbi"), emit: tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    /usr/local/bin/run_clair3.sh \
    --bam_fn=${input} \
    --ref_fn=${fasta} \
    --threads=${task.cpus} \
    --platform=${params.platform} \
    --model_path="/usr/local/bin/models/${params.model_name}" \
    --output="${meta.id}" \
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /usr/local/bin/run_clair3.sh --version | sed 's/Clair3 //' )
    END_VERSIONS
    """
}
