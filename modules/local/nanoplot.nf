process NANOPLOT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::nanoplot=1.38.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.38.0--pyhdfd78af_0' :
        'quay.io/biocontainers/nanoplot:1.38.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("$output_html"), emit: html
    //tuple val(meta), path("$output_png") , emit: png
    tuple val(meta), path("$output_txt") , emit: txt
    tuple val(meta), path("$output_log") , emit: log
    path  "versions.yml"           , emit: versions

    script:
    //def args = task.ext.args ?: ''
    // $options.args \\
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
        -t $task.cpus \\
        $input_file \\
        -o $output_dir
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
