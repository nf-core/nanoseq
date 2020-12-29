// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
def options          = initOptions(params.options)

process FASTQC {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "quay.io/biocontainers/fastqc:0.11.9--0"
    
    input:
    tuple val(sample), path(fastq)
    
    output:    
    path "*.zip"                , emit: zip
    path "*.html"               , emit: html
    path "*.version.txt"        , emit: version

    script:
    """
    [ ! -f  ${sample}.fastq.gz ] && ln -s $fastq ${sample}.fastq.gz
    fastqc \\
        -q \\
        -t $task.cpus \\
        ${sample}.fastq.gz
    fastqc --version > fastqc.version.txt
    """
}
