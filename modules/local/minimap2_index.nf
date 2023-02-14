process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::minimap2=2.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3' :
        'quay.io/biocontainers/minimap2:2.17--hed695b0_3' }"

    input:
    tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), val(annotation_str)

    output:
    tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path("*.mmi"), val(annotation_str), emit: index
    path "versions.yml"                                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "-ax map-ont" : "-ax splice"
    def kmer      = (params.protocol == 'directRNA') ? "-k14" : ""
    def stranded  = (params.stranded || params.protocol == 'directRNA') ? "-uf" : ""
    def junctions = (params.protocol != 'DNA' && bed) ? "--junc-bed ${file(bed)}" : ""
    """
    minimap2 \\
        $preset \\
        $kmer \\
        $stranded \\
        $junctions \\
        -t $task.cpus \\
        -d ${fasta}.mmi \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
