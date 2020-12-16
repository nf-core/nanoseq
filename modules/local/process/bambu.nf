// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BAMBU {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:sample) }

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=1.0.0" : null)
    container "docker.io/bioconda::bioconductor-bambu=1.0.0" //need a multitool container for bambu and r-base on quay hub 

    input:
    tuple path(fasta), path(gtf)
    path bams from ch_sortbam_quant.collect{ it[1] }
    
    output:
    path "counts_gene.txt"          , emit: gene_counts
    path "counts_transcript.txt"    , emit: transcript_counts
    path "extended_annotations.gtf" , emit: extended_gtf

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=$gtf \\
        --fasta=$fasta $bams
    """
}
