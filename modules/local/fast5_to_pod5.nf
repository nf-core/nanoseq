process FAST5_TO_POD5 {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.0.8 bioconda::bioconductor-bsgenome=1.66.0"
    container "docker.io/yuukiiwa/pod5:0.2.4"

    input:
    tuple val(meta), path(input_path)

    output:
    tuple val(meta), path("pod5/")    , emit: pod5

    when:
    task.ext.when == null || task.ext.when

    script:
    output_name = "pod5/converted.pod5"
    """
    pod5 convert fast5 $input_path --output $output_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$(echo \$(pod5 --version 2>&1) | sed -r 's/..............//')
    END_VERSIONS
    """
}
