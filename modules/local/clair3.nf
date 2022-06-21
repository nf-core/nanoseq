process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:0.1.10--hdfd78af_0' :
        'quay.io/biocontainers/clair3:0.1.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(input), path(index)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}/*")        , emit: vcf
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    /usr/local/bin/run_clair3.sh \
    --bam_fn=$input \
    --ref_fn=$fasta \
    --threads=$task.cpus \
    --platform="ont" \
    --model_path="/usr/local/bin/models/${params.clair3_model}" \
    --output="${meta.id}" \
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /usr/local/bin/run_clair3.sh --version | sed 's/ /,/' )
    END_VERSIONS
    """
}
