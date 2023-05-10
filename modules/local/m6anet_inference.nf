process M6ANET_INFERENCE {
    tag "$meta.id"
    echo true
    label 'process_medium'

//  conda     (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null) // need to get xpore onto conda
    container "docker.io/yuukiiwa/m6anet:1.0"

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
    m6anet-run_inference --input_dir $input_dir --out_dir $out_dir  --batch_size 512 --n_processes $task.cpus --num_iterations 5 --device cpu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        m6anet: \$( echo 'm6anet 1.0' )
    END_VERSIONS
    """
}
