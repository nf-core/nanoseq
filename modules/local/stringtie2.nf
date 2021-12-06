// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process STRINGTIE2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.1.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/stringtie:2.1.4--h7e0af3c_0"
    } else {
        container "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"
    }

    input:
    tuple val(meta), path(fasta), path(gtf), path(bam)

    output:
    path "*.stringtie.gtf"       , emit: stringtie_gtf
    path  "versions.yml"         , emit: versions

    script:
    """
    stringtie \\
        -L \\
        -G $gtf \\
        -o ${meta.id}.stringtie.gtf $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
