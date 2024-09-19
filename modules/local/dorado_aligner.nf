process DORADO_ALIGNER {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/ontresearch/dorado"

    input:
    tuple val(meta), path(mod_bam)
    path fasta

    output:
    tuple val(meta), path("aligned_sorted.bam"), path("*.bai")  , emit: aligned_bam
    path "versions.yml"                                  , emit: versions

    script:
    """
    dorado aligner --mm2-preset map-ont $fasta $mod_bam > aligned.bam && samtools sort aligned.bam -o aligned_sorted.bam && samtools index aligned_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS
    """
}

