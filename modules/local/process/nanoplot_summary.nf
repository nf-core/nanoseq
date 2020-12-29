// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOPLOT_SUMMARY {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    conda     (params.enable_conda ? "bioconda::nanoplot=1.32.1" : null)
    container "quay.io/biocontainers/nanoplot:1.30.1--py_0"

    when:
    !params.skip_basecalling && !params.skip_qc && !params.skip_nanoplot

    input:
    path summary_txt
    
    output:
    path "*.png"         , emit: png
    path "*.html"        , emit: html 
    path "*.txt"         , emit: txt  
    path "*.log"         , emit: log

    script:
    """
    NanoPlot -t $task.cpus --summary $summary_txt
    """
}
