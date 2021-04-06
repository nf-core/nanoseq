// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

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
    path "dexseq.version.txt"   , emit: dexseq_version
    path "drimseq.version.txt"  , emit: drimseq_version
    path "stager.version.txt"   , emit: stager_version

    script:
    def software    = getSoftwareName(task.process)
    def label_lower = params.multiqc_label.toLowerCase()
    def label_upper = params.multiqc_label.toUpperCase()
    """
    run_dexseq.r $params.quantification_method $counts
    Rscript -e "library(DEXSeq); write(x=as.character(packageVersion('DEXSeq')), file='dexseq.version.txt')"
    Rscript -e "library(DRIMSeq); write(x=as.character(packageVersion('DRIMSeq')), file='drimseq.version.txt')"
    Rscript -e "library(stageR); write(x=as.character(packageVersion('stageR')), file='stager.version.txt')"
    """
}
