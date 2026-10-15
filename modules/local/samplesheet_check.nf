process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    path samplesheet
    path input_path   // Directory used to resolve relative 'input_file' entries, or [] if not needed

    output:
    path '*.csv', emit: csv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    script: // This script is bundled with the pipeline, in nf-core/nanoseq/bin/
    def updated_path = input_path ? "\$(realpath ${input_path})" : 'not_changed'
    """
    check_samplesheet.py \\
        $samplesheet \\
        $updated_path \\
        samplesheet.valid.csv
    """
}
