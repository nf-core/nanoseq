// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options       = [:]
params.multiqc_label = ''
def options          = initOptions(params.options)

process DESEQ2 {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.3 bioconda::bioconductor-deseq2=1.28.0" : null)
    container "quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
    
    input:
    path counts
    
    output:    
    path "*.txt"                , emit: deseq2_txt
    path "versions.yml"         , emit: versions

    script:
    """
    run_deseq2.r $params.quantification_method $counts

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
