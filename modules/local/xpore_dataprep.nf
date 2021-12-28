process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

//  conda     (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null) // need to get xpore onto conda
    container "docker.io/yuukiiwa/xpore:2.1"

    input:
    tuple val(meta), path(genome), path(gtf), path(eventalign), path(nanopolish_summary)

    output:
    tuple val(meta), path("$meta.id"), emit: dataprep_outputs

    script:
    """
    xpore dataprep \\
    --eventalign $eventalign \\
    --out_dir $meta.id \\
    --n_processes $task.cpus \\
    --genome --gtf_or_gff $gtf --transcript_fasta $genome
    """
}
