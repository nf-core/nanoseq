process XPORE_DIFFMOD {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    // conditions[i] is the condition label for dataprep_dirs[i]
    tuple val(conditions), path(dataprep_dirs)

    output:
    path "diffmod*", emit: diffmod_outputs
    tuple val("${task.process}"), val('xpore'), eval("xpore --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_xpore, topic: versions

    script:
    def args    = task.ext.args ?: ''
    def dirs    = dataprep_dirs instanceof List ? dataprep_dirs : [dataprep_dirs]
    def samples = [dirs, conditions].transpose().collect { pair -> "\$PWD/${pair[0]},${pair[1]}" }.join(' ')
    """
    create_yml.py diffmod_config.yml $samples
    xpore diffmod \\
        --config diffmod_config.yml \\
        --n_processes $task.cpus \\
        $args
    """

    stub:
    """
    touch diffmod_config.yml
    mkdir -p diffmod_outputs
    touch diffmod_outputs/diffmod.table
    touch diffmod_outputs/diffmod.log
    """
}
