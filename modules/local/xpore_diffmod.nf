process XPORE_DIFFMOD {
    label 'process_medium'

//  conda     (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null) // need to get xpore onto conda
    container "docker.io/yuukiiwa/xpore:2.1"

    input:
    val dataprep_dirs

    output:
    path "diffmod*", emit: diffmod_outputs
    path "versions.yml"        , emit: versions

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
