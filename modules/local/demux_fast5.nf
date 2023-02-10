process DEMUX_FAST5 {
    label 'process_medium'

    conda "bioconda::ont-fast5-api=4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-fast5-api:4.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/ont-fast5-api:4.0.0--pyhdfd78af_0' }"

    input:
    path(input_path), stageAs: 'input_path/*'
    tuple val(meta), path(input_summary)

    output:
    path "demultiplexed_fast5/*"   , emit: fast5
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fast5_dir_path = workflow.profile.contains('test') ? "input_path" : "$input_path"
    """
    demux_fast5 \\
    --input  $fast5_dir_path \\
    --save_path ./demultiplexed_fast5 \\
    --summary_file $input_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        demux_fast5: \$(python -c 'import ont_fast5_api;print(ont_fast5_api.__version__)')
    END_VERSIONS
    """
}
