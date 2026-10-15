process DEXSEQ {
    label "process_medium"

    conda "conda-forge::r-base=4.0.2 bioconda::bioconductor-dexseq=1.36.0 bioconda::bioconductor-drimseq=1.18.0 bioconda::bioconductor-stager=1.12.0"
    // TODO: replace with a multi-tool (mulled/Seqera Containers) image for r-base, dexseq, stager and drimseq
    container "docker.io/yuukiiwa/nanoseq:dexseq"

    input:
    path counts
    val quantification_method

    output:
    path "*.txt", emit: dexseq_txt
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | sed '1!d;s/^R version //;s/ .*//'"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('bioconductor-dexseq'), eval("Rscript -e \"cat(as.character(packageVersion('DEXSeq')))\""), emit: versions_dexseq, topic: versions
    tuple val("${task.process}"), val('bioconductor-drimseq'), eval("Rscript -e \"cat(as.character(packageVersion('DRIMSeq')))\""), emit: versions_drimseq, topic: versions
    tuple val("${task.process}"), val('bioconductor-stager'), eval("Rscript -e \"cat(as.character(packageVersion('stageR')))\""), emit: versions_stager, topic: versions

    script:
    """
    run_dexseq.r $quantification_method $counts
    """

    stub:
    """
    touch dexseq.results.txt
    """
}
