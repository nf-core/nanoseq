// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
def options          = initOptions(params.options)

process FASTQC {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "quay.io/biocontainers/fastqc:0.11.9--0"
    
    input:
    tuple val(meta), path(fastq)
    
    output:    
    path "*.zip"                , emit: zip
    path "*.html"               , emit: html
    path "*.version.txt"        , emit: version

    script:
    """
    [ ! -f  ${meta.id}.fastq.gz ] && ln -s $fastq ${meta.id}.fastq.gz
    fastqc \\
        -q \\
        -t $task.cpus \\
        ${meta.id}.fastq.gz
    fastqc --version > fastqc.version.txt
    """
}
