// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]

/*
 * Convert GTF file to BED format
 */
process GTF2BED {
    label 'process_low'
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::perl=5.26.2" : null)
    container "quay.io/biocontainers/perl:5.26.2"

    input:
    tuple path(gtf), val(name)

    output:
    tuple path('*.bed'), val(name) , emit: gtf_bed
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
