// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NANOPLOT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::nanoplot=1.38.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanoplot:1.38.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("$output_html"), emit: html
    tuple val(meta), path("$output_png") , emit: png
    tuple val(meta), path("$output_txt") , emit: txt
    tuple val(meta), path("$output_log") , emit: log
    path  "versions.yml"           , emit: versions

    script:
    def input_file = ("$ontfile".endsWith(".fastq.gz")) ? "--fastq ${ontfile}" :
                     ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    def output_dir = ("$ontfile".endsWith(".fastq.gz")) ? "fastq/${meta.id}" :
                     ("$ontfile".endsWith(".txt")) ? "summary" : ''
    output_html = output_dir+"/*.html"
    output_png  = output_dir+"/*.png"
    output_txt  = output_dir+"/*.txt"
    output_log  = output_dir+"/*.log"
    """
    NanoPlot \\
        $options.args \\
        -t $task.cpus \\
        $input_file \\
        -o $output_dir
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
