// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

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
    container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    if ({ task.exitStatus in [255] }) { errorStrategy = 'ignore' }

    input:
    tuple val(meta), path(sizes), path(bedgraph)
    
    output:
    tuple val(meta), path(sizes), path("*.bigWig"), emit: bigwig
    path "*.version.txt"             , emit: version

    script:
    """
    bedGraphToBigWig $bedgraph $sizes ${meta.id}.bigWig
    echo $VERSION > bedGraphToBigWig.version.txt
    """
}
