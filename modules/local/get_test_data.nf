process GET_TEST_DATA {
    label "process_single"

    container "docker.io/yuukiiwa/git:latest"

    output:
    path "test-datasets/fast5/$barcoded/*"        , emit: ch_input_fast5s_path
    path "test-datasets/modification_fast5_fastq/", emit: ch_input_dir_path
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    barcoded = (workflow.profile.contains('test_bc_nodx') || workflow.profile.contains('rnamod')) ? "nonbarcoded" : "barcoded"
    """
    git clone https://github.com/nf-core/test-datasets.git --branch nanoseq --single-branch

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        git: \$(echo \$(git --version | sed 's/git version //; s/ .*\$//')
    """
}
