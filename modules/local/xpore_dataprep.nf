process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "bioconda::xpore=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh5e36f6f_0' :
        'quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(eventalign), path(nanopolish_summary)

    output:
    tuple val(meta), path("$meta.id"), emit: dataprep_outputs
    path "versions.yml"        , emit: versions

    script:
    """
    xpore dataprep \\
    --eventalign $eventalign \\
    --out_dir $meta.id \\
    --n_processes $task.cpus \\
    --genome --gtf_or_gff $gtf --transcript_fasta $genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$( xpore --version | sed -e 's/xpore version //g' )
    END_VERSIONS
    """
}
