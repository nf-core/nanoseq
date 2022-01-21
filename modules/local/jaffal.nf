process JAFFAL {
    tag "$meta.id"
    echo true
    label 'process_medium'

    conda     (params.enable_conda ? "bioconda::jaffa=2.0.0" : null)
    container "docker.io/yuukiiwa/jaffa:2.0"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/jaffa:2.00--hdfd78af_1' :
    //    'quay.io/biocontainers/jaffa:2.00--hdfd78af_1' }"//tried three biocontainers, all of them got command not found for minimap2 

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
    echo 'jaffa 2.0' > jaffal_version.txt
    """
}
