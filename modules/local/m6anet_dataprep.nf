process M6ANET_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::m6anet=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/m6anet:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/m6anet:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(eventalign), path(nanopolish_summary)

    output:
    tuple val(meta), path("${prefix}"), emit: dataprep_outputs
    tuple val("${task.process}"), val('m6anet'), eval("m6anet --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_m6anet, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    m6anet dataprep \\
        --eventalign $eventalign \\
        --out_dir $prefix \\
        --n_processes $task.cpus
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    touch $prefix/data.json
    touch $prefix/data.info
    touch $prefix/data.log
    """
}
