process DORADO {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/ontresearch/dorado"

    input:
    tuple val(meta), path(pod5_path)
    val dorado_device
    val dorado_model

    output:
    tuple val(meta), path("*.fastq.gz")  , emit: fastq
    path "versions.yml"                  , emit: versions

    script:
    """
    dorado download --model $dorado_model
    dorado basecaller $dorado_model $pod5_path --device $dorado_device --emit-fastq > basecall.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS

    gzip basecall.fastq
    """
}

