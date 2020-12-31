// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BAMBU {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=1.0.2 bioconda::bioconductor-bsgenome=1.58.0" : null)
    container "docker.io/yuukiiwa/nanoseq:bambu_bsgenome"

    input:
    tuple path(fasta), path(gtf)
    path bams
    
    output:
    path "counts_gene.txt"          , emit: ch_gene_counts
    path "counts_transcript.txt"    , emit: ch_transcript_counts
    path "extended_annotations.gtf" , emit: extended_gtf
    path "bambu.version.txt"        , emit: bambu_version
    path "r.version.txt"            , emit: r_version

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams
    Rscript -e "library(bambu); write(x=as.character(packageVersion('bambu')), file='bambu.version.txt')"
    echo \$(R --version 2>&1) > r.version.txt
    """
}
