process DORADO_BASECALLER {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/ontresearch/dorado"

    input:
    tuple val(meta), path(pod5_path)
    val dorado_device
    val dorado_model

    output:
    tuple val(meta), path("basecall*")  , emit: dorado_out
    path "versions.yml"                  , emit: versions

    script:
    def emit_args = (params.dorado_modification == null) ? " --emit-fastq > basecall.fastq && gzip basecall.fastq" : " --modified-bases $params.dorado_modification > basecall.bam"
    """
    dorado download --model $dorado_model
    dorado basecaller $dorado_model $pod5_path --device $dorado_device $emit_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}

