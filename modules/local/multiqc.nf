// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"

    input:
    path ch_multiqc_config
    path ch_multiqc_custom_config
    path ch_pycoqc_multiqc
    path ch_fastqc_multiqc
    path ch_sortbam_stats_multiqc
    path ch_featurecounts_gene_multiqc
    path ch_featurecounts_transcript_multiqc
    path software_versions_yaml
    path workflow_summary
    
    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f . $custom_config

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

