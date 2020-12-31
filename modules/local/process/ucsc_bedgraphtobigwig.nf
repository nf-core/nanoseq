// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '377'

process UCSC_BEDGRAPHTOBIGWIG {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::ucsc-bedgraphtobigwig=377" : null)
    container "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1"
    if ({ task.exitStatus in [255] }) { errorStrategy = 'ignore' }

    when:
    !params.skip_alignment && !params.skip_bigwig && (params.protocol == 'directRNA' || params.protocol == 'cDNA')

    input:
    tuple val(sample), path(sizes), path(bedgraph)
    
    output:
    tuple val(sample), path(sizes), path("*.bigWig"), emit: bigwig

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $sizes \\
        ${sample}.bigWig
    echo $VERSION > bedGraphToBigWig.version.txt
    """
}
