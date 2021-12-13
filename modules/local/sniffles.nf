// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNIFFLES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::sniffles=1.0.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/sniffles:1.0.12--h8b12597_1'
    } else {
        container 'quay.io/biocontainers/sniffles:1.0.12--h8b12597_1'
    }

    input:
    tuple val(meta), path(sizes), val(is_transcripts), path(bam), path(bai)


    output:
    path "*_sniffles.vcf"               , emit: sv_calls // vcf files
    path "versions.yml"                 , emit: versions


    script:
    """
    sniffles \
        -m  $bam \
        -v ${meta.id}_sniffles.vcf \
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}