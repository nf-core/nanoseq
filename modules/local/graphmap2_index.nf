process GRAPHMAP2_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::graphmap=0.6.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'quay.io/biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), val(annotation_str)

    output:
    tuple path(fasta), path(sizes), path(gtf), val(bed), val(is_transcripts), path("*.gmidx"), val(annotation_str), emit: index
    path "versions.yml"                                                                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
    def junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf $gtf" : ""
    """
    graphmap2 \\
        align \\
        $preset \\
        $junctions \\
        -t $task.cpus \\
        -I \\
        -r $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphmap2: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
