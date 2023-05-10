process GRAPHMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::graphmap=0.6.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'quay.io/biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    tuple val(meta), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)

    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
    def junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf $gtf" : ""
    """
    graphmap2 \\
        align \\
        $preset \\
        $junctions \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $fastq \\
        -o ${meta.id}.sam \\
        --extcigar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphmap2: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
