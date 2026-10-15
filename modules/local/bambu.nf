process BAMBU {
    label 'process_medium'

    conda "conda-forge::r-base=4.3.1 bioconda::bioconductor-bambu=3.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.4.0--r43hf17093f_0' :
        'quay.io/biocontainers/bioconductor-bambu:3.4.0--r43hf17093f_0' }"

    input:
    tuple path(fasta), path(gtf)
    path bams

    output:
    path "counts_gene.txt"         , emit: gene_counts
    path "counts_transcript.txt"   , emit: transcript_counts
    path "extended_annotations.gtf", emit: extended_gtf
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | sed '1!d;s/^R version //;s/ .*//'"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('bioconductor-bambu'), eval("Rscript -e \"cat(as.character(packageVersion('bambu')))\""), emit: versions_bambu, topic: versions
    tuple val("${task.process}"), val('bioconductor-bsgenome'), eval("Rscript -e \"cat(as.character(packageVersion('BSgenome')))\""), emit: versions_bsgenome, topic: versions

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams
    """
}
