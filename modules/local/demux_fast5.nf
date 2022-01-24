process DEMUX_FAST5 {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    conda     (params.enable_conda ? "bioconda:ont-fast5-api:4.0.0--pyhdfd78af_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ont-fast5-api:4.0.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/ont-fast5-api:4.0.0--pyhdfd78af_0"
    }

    input:
    path(input_path), stageAs: 'input_path/*'
    tuple val(meta), path(input_summary)

    output:
    path "demultiplexed_fast5/*"   , emit: fast5
    path "versions.yml"            , emit: versions

    script:
    """
    demux_fast5 \\
    --input  input_path \\
    --save_path ./demultiplexed_fast5 \\
    --summary_file $input_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    demux_fast5: \$(echo \$(python -c\'import ont_fast5_api;print(ont_fast5_api.__version__)\'))
    END_VERSIONS
    """
}
