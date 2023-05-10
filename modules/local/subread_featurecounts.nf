process SUBREAD_FEATURECOUNTS {
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "bioconda::subread=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'quay.io/biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    path gtf
    path bams

    output:
    path "counts_gene.txt"              , emit: gene_counts
    path "counts_transcript.txt"        , emit: transcript_counts
    path "counts_gene.txt.summary"      , emit: featurecounts_gene_multiqc
    path "counts_transcript.txt.summary", emit: featurecounts_transcript_multiqc
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
    "${task.process}":
        featureCounts: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}
