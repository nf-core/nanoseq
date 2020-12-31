// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BAM_RENAME {
    tag "$sample"

    input:
    tuple val(sample), path(bam)
    
    output:
    tuple val(sample), path("*.bam"), emit: sortbam_quant

    script:
    """
    [ ! -f ${sample}.bam ] && ln -s $bam ${sample}.bam
    """
}
