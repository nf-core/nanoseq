process DEXSEQ {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1c0ca9d50fbf5c4a72419bca15a95fe8f3bbf741f5b203d03a100cf384d5172/data' :
        'community.wave.seqera.io/library/r-base_bioconductor-dexseq_bioconductor-drimseq_bioconductor-stager:ad53f230703306c9' }"

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
