// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process STRINGTIE2 {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:sample) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.4" : null)
    container "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"

    input:
    tuple val(sample), path(fasta), path(gtf), path(bam)
    
    output:
    path "*.stringtie.gtf"   , emit: stringtie_gtf
    path  "*.version.txt"    , emit: version

    script:
    """
    stringtie \\
        -L \\
        -G $gtf \\
        -o ${sample}.stringtie.gtf $bam
    stringtie --version > ${software}.version.txt
    """
}
