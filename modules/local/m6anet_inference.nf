process M6ANET_INFERENCE {
    tag "$meta.id"
    label 'process_medium'
    debug true

    conda "bioconda::m6anet=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/m6anet:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/m6anet:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_dir)

    output:
    path "${prefix}_results", emit: m6anet_outputs
    tuple val("${task.process}"), val('m6anet'), eval("m6anet --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_m6anet, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    m6anet inference \\
        --input_dir $input_dir \\
        --out_dir ${prefix}_results \\
        --batch_size 512 \\
        --n_processes $task.cpus \\
        --num_iterations 5 \\
        --device cpu
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_results
    touch ${prefix}_results/data.site_proba.csv
    touch ${prefix}_results/data.indiv_proba.csv
    """
}
