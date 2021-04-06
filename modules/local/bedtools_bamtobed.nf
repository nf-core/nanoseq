// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// params.options = [:]
// def options    = initOptions(params.options)

process BEDTOOLS_BAMBED {
    label 'process_medium'

    conda     (params.enable_conda ? "bioconda::bedtools=2.29.2" : null)
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"

    when:
    !params.skip_alignment && !params.skip_bigbed && (params.protocol == 'directRNA' || params.protocol == 'cDNA')

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(bam), path(bai)
    
    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.bed12"), emit: bed12
    path "*.version.txt"                                                , emit: version

    script:
    """
    bedtools \\
        bamtobed \\
        -bed12 \\
        -cigar \\
        -i ${bam[0]} \\
        | bedtools sort > ${meta.id}.bed12
    bedtools --version | sed -e "s/bedtools v//g" > bedtools.version.txt
    """
}
