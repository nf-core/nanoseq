process M6ANET_INFERENCE {
    tag "$meta.id"
    echo true
    label 'process_medium'

//  conda     (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null) // need to get xpore onto conda
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/m6anet:2.0.2--pyhdfd78af_0' :
        'quay.io/biocontainers/m6anet:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_dir)

    output:
    path "*"           , emit: m6anet_outputs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_dir = meta.id+"_results"
    """
    m6anet inference --input_dir $input_dir --out_dir $out_dir  --batch_size 512 --n_processes $task.cpus --num_iterations 5 --device cpu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        m6anet: \$( echo 'm6anet 2.0.2' )
    END_VERSIONS
    """
}
