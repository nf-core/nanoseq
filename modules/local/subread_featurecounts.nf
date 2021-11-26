// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SUBREAD_FEATURECOUNTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::subread=2.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0"
    } else {
        container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    }

    input:
    path gtf
    path bams

    output:
    path "counts_gene.txt"               , emit: gene_counts
    path "counts_transcript.txt"         , emit: transcript_counts
    path "counts_gene.txt.summary"       , emit: featurecounts_gene_multiqc
    path "counts_transcript.txt.summary" , emit: featurecounts_transcript_multiqc
    path "versions.yml"                  , emit: versions

    script:
    """
    featureCounts \\
        -L \\
        -O \\
        -f \\
        -g gene_id \\
        -t exon \\
        -T $task.cpus \\
        -a $gtf \\
        -o counts_gene.txt \\
        $bams

    featureCounts \\
        -L \\
        -O \\
        -f \\
        --primary \\
        --fraction \\
        -F GTF \\
        -g transcript_id \\
        -t transcript \\
        --extraAttributes gene_id \\
        -T $task.cpus \\
        -a $gtf \\
        -o counts_transcript.txt \\
        $bams

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}
