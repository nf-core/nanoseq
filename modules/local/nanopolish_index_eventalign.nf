process NANOPOLISH_INDEX_EVENTALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::nanopolish==0.13.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.13.2--he3b7ca5_2' :
        'quay.io/biocontainers/nanopolish:0.13.2--he3b7ca5_2' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fast5), path(fastq), path(bam), path(bai)

    output:
    tuple val(meta), path(genome), path(gtf), path("*eventalign.txt"), path("*summary.txt"), emit: nanopolish_outputs
    path "versions.yml"        , emit: versions

    script:
    sample_summary = "$meta.id" +"_summary.txt"
    sample_eventalign = "$meta.id" +"_eventalign.txt"
    """
    nanopolish index -d $fast5 $fastq
    nanopolish eventalign  --reads $fastq --bam $bam --genome $genome --scale-events --signal-index --summary $sample_summary --threads $task.cpus > $sample_eventalign

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$( nanopolish --version | sed -e 's/nanopolish version //g' )
    END_VERSIONS
    """
}
