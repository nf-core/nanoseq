process JAFFAL {
    tag "$meta.id"
    label 'process_medium'
    debug true

    conda "bioconda::jaffa=2.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jaffa:2.5--hdfd78af_0' :
        'quay.io/biocontainers/jaffa:2.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path jaffal_ref_dir

    output:
    tuple val(meta), path("*.fasta"), emit: jaffal_fastq
    path "*.csv"                    , emit: jaffal_results
    // JAFFA has no command-line version option; keep in sync with the container tag
    tuple val("${task.process}"), val('jaffa'), val('2.5'), emit: versions_jaffa, topic: versions

    script:
    """
    bpipe run -p refBase=$jaffal_ref_dir $jaffal_ref_dir/JAFFAL.groovy $fastq
    """
}
