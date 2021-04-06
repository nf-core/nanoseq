// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PYCOQC {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    container "quay.io/biocontainers/pycoqc:2.5.0.21--py_0"

    when:
    !params.skip_basecalling && !params.skip_qc && !params.skip_pycoqc

    input:
    tuple val(meta), path(summary_txt)
    
    output:
    path "*.html"                      , emit: html
    path "*.json"                      , emit: pycoqc_multiqc
    path "*.version.txt"               , emit: version

    script:
    """
    pycoQC -f $summary_txt -o pycoqc.html -j pycoqc.json
    pycoQC --version &> pycoqc.version.txt
    """
}
