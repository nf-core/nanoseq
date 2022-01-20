process BAMBU {
    label 'process_medium'

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=2.0.0 bioconda::bioconductor-bsgenome=1.62.0" : null)
    container "docker.io/yuukiiwa/nanoseq:bambu_bsgenome" //not on biocontainers; does not have a singularity container

    input:
    tuple path(fasta), path(gtf)
    path bams

    output:
    path "counts_gene.txt"          , emit: ch_gene_counts
    path "counts_transcript.txt"    , emit: ch_transcript_counts
    path "extended_annotations.gtf" , emit: extended_gtf
    path "versions.yml"             , emit: versions

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}
