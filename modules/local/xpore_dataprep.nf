process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::xpore=2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(eventalign), path(nanopolish_summary)

    output:
    tuple val(meta), path("${prefix}"), emit: dataprep_outputs
    tuple val("${task.process}"), val('xpore'), eval("xpore --version 2>&1 | sed '1!d;s/.* //'"), emit: versions_xpore, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    xpore dataprep \\
        --eventalign $eventalign \\
        --out_dir $prefix \\
        --n_processes $task.cpus \\
        --genome \\
        --gtf_or_gff $gtf \\
        --transcript_fasta $genome
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    touch $prefix/data.json
    touch $prefix/data.index
    touch $prefix/data.readcount
    touch $prefix/data.log
    """
}
