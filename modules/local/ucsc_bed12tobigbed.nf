// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '377'

process UCSC_BED12TOBIGBED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'bigBed', publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::ucsc-bedtobigbed=377" : null)
    container "quay.io/biocontainers/ucsc-bedtobigbed:377--h446ed27_1"
    if ({ task.exitStatus in [255] }) { errorStrategy = 'ignore' }

    when:
    !params.skip_alignment && !params.skip_bigbed && (params.protocol == 'directRNA' || params.protocol == 'cDNA')

    input:
    tuple val(meta), path(sizes),  path(bed12)

    output:
    tuple val(meta), path(sizes),  path("*.bigBed"), emit: bigbed
    path "versions.yml"                            , emit: versions

    script:
    """
    bedToBigBed \\
        $bed12 \\
        $sizes \\
        ${meta.id}.bigBed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
