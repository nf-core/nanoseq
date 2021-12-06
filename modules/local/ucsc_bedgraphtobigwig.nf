// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
// def options    = initOptions(params.options)  //this line gives "Cannot get property 'args' on null object"

def VERSION = '377'

process UCSC_BEDGRAPHTOBIGWIG {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:377--h446ed27_1"
    } else {
        container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    }

    input:
    tuple val(meta), path(sizes), path(bedgraph)

    output:
    tuple val(meta), path(sizes), path("*.bigWig"), emit: bigwig
    path "versions.yml"                           , emit: versions

    script:
    """
    bedGraphToBigWig $bedgraph $sizes ${meta.id}.bigWig

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
