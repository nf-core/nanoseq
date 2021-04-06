// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SUBREAD_FEATURECOUNTS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::subread=2.0.1" : null)
    container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"

    input:
    path gtf
    path bams
    
    output:
    path "counts_gene.txt"               , emit: gene_counts
    path "counts_transcript.txt"         , emit: transcript_counts
    path "counts_gene.txt.summary"       , emit: featurecounts_gene_multiqc
    path "counts_transcript.txt.summary" , emit: featurecounts_transcript_multiqc
    path "*.version.txt"                 , emit: version

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
    echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" > featureCounts.version.txt
    """
}
