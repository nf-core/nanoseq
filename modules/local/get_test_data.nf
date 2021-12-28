process GET_TEST_DATA {
    container "docker.io/yuukiiwa/git:latest"

    output:    
    path "test-datasets/$subdir/$barcoded"  , emit: ch_input_path

    script:
    subdir = workflow.profile.contains('test_nobc_nodx_rnamod') ? "modification_fast5_fastq" : "fast5"
    barcoded = (workflow.profile.contains('test_bc_nodx') || workflow.profile.contains('test_nobc_nodx_rnamod')) ? "nonbarcoded" : "barcoded"
    """
    git clone https://github.com/nf-core/test-datasets.git --branch nanoseq --single-branch
    """
}
