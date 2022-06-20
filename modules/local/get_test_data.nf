process GET_TEST_DATA {
    container "docker.io/yuukiiwa/git:latest"

    output:
    path "test-datasets/fast5/$barcoded/*"  , emit: ch_input_fast5s_path
    path "test-datasets/modification_fast5_fastq/"   , emit: ch_input_dir_path

    script:
    barcoded = (workflow.profile.contains('test_bc_nodx') || workflow.profile.contains('rnamod')) ? "nonbarcoded" : "barcoded"
    """
    git clone https://github.com/nf-core/test-datasets.git --branch nanoseq --single-branch
    """
}
