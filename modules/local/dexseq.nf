// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options       = [:]
params.multiqc_label = ''
def options          = initOptions(params.options)

process DEXSEQ {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0" : null)
    container "docker.io/yuukiiwa/nanoseq:dexseq" 
    // need a multitool container for r-base, dexseq, stager, drimseq and on quay hub

    input:
    path counts
    
    output:    
    path "*.txt"                , emit: dexseq_txt
    path "versions.yml"         , emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def label_lower = params.multiqc_label.toLowerCase()
    def label_upper = params.multiqc_label.toUpperCase()
    """
    run_dexseq.r $params.quantification_method $counts

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DEXSeq); cat(as.character(packageVersion('DEXSeq')))")
        bioconductor-deseq2: \$(Rscript -e "library(DRIMSeq); cat(as.character(packageVersion('DRIMSeq')))")
        bioconductor-deseq2: \$(Rscript -e "library(stageR); cat(as.character(packageVersion('stageR')))")
    END_VERSIONS
    """
}
