// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3"
    } else {
        container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
    }

    input:
    tuple val(meta), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)

    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam
    path "versions.yml" , emit: versions

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
        $index \\
        $fastq > ${meta.id}.sam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
