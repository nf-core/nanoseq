// Import generic module functions
include { saveFiles;; getProcessName } from './functions'

process BAM_RENAME {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.bam"), emit: sortbam_quant

    script:
    """
    [ ! -f ${meta.id}.bam ] && ln -s $bam ${meta.id}.bam
    """
}
