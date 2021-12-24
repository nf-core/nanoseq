process JAFFAL {
    tag "$meta.id"
    echo true
    label 'process_medium'

//    conda     (params.enable_conda ? "bioconda::nanolyse=1.2.0" : null)   //need conda container
    container "docker.io/yuukiiwa/jaffa:2.0"

    input:
    tuple val(meta), path(fastq)
    path jaffal_ref_dir

    output:
    tuple val(meta), path("*.fasta") ,emit: jaffal_fastq
    path "*.csv"                     ,emit: jaffal_results
    path "*_version.txt"             ,emit: version

    script:
    """
    bpipe run -p refBase=$jaffal_ref_dir $jaffal_ref_dir/JAFFAL.groovy $fastq
    echo '2.0' > jaffal_version.txt
    """
}
