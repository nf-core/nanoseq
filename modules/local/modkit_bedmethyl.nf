process DORADO {
    tag "$meta.id"
    label 'process_medium'

    container "modkit"

    input:
    tuple val(meta), path(mod_bam)

    output:
    tuple val(meta), path("*.bed") , emit: bedmethyl
    path "versions.yml"                  , emit: versions

    script:
    bedmethyl = "$meta.id" +".bed"
    """
    modkit pileup $mod_bam $bedmethyl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS

    gzip basecall.fastq
    """
}

