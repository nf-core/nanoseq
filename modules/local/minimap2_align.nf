process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::minimap2=2.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3' :
        'quay.io/biocontainers/minimap2:2.17--hed695b0_3' }"

    input:
    tuple val(meta), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)

    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "-ax map-ont" : "-ax splice"
    def kmer      = (params.protocol == 'directRNA') ? "-k14" : ""
    def stranded  = (params.stranded || params.protocol == 'directRNA') ? "-uf" : ""
    def junctions = (params.protocol != 'DNA' && bed) ? "--junc-bed ${file(bed)}" : ""
    def md        = (params.call_variants && params.protocol == 'DNA') ? "--MD" : ""
    """
    minimap2 \\
        $preset \\
        $kmer \\
        $stranded \\
        $junctions \\
        $md \\
        -t $task.cpus \\
        $index \\
        $fastq > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
