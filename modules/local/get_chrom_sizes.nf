// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

/*
 * Get chromosome sizes from a fasta file
 */
process GET_CHROM_SIZES {
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, publish_dir:"genome") }

    conda     (params.enable_conda ? "bioconda::samtools=1.10" : null)
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    input:
    tuple path(fasta), val(name)

    output:
    tuple path('*.sizes'), val(name) , emit: sizes
    //  path '*.fai'                     , emit: fai
    path  "versions.yml"            , emit: versions

    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
