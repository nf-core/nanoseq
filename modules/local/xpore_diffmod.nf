process XPORE_DIFFMOD {
    label 'process_medium'

    conda "bioconda::xpore=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    val dataprep_dirs

    output:
    path "diffmod*"    , emit: diffmod_outputs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    diffmod_config = "--config $workflow.workDir/*/*/diffmod_config.yml"
    """
    create_yml.py diffmod_config.yml $dataprep_dirs
    xpore diffmod $diffmod_config --n_processes $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$( xpore --version | sed -e 's/xpore version //g' )
    END_VERSIONS
    """
}
