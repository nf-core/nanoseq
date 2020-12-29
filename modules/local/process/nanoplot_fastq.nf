// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOPLOT_FASTQ {
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::nanoplot=1.32.1" : null)
    container "quay.io/biocontainers/nanoplot:1.30.1--py_0"

    when:
    !params.skip_qc && !params.skip_nanoplot

    input:
    tuple val(sample), path(fastq)
    
    output:
    path "*.png"               , emit: png
    path "*.html"              , emit: html
    path "*.txt"               , emit: txt
    path "*.log"               , emit: log
    path "*.version.txt"       , emit: version

    script:
    """
    NanoPlot -t $task.cpus --fastq $fastq
    NanoPlot --version &> nanoplot.version.txt
    """
}
