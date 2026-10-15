process DESEQ2 {
    label "process_medium"

    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-deseq2=1.28.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path counts
    val quantification_method

    output:
    path "*.txt", emit: deseq2_txt
    tuple val("${task.process}"), val('r-base'), eval("R --version 2>&1 | sed '1!d;s/^R version //;s/ .*//'"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('bioconductor-deseq2'), eval("Rscript -e \"cat(as.character(packageVersion('DESeq2')))\""), emit: versions_deseq2, topic: versions

    script:
    """
    run_deseq2.r $quantification_method $counts
    """

    stub:
    """
    touch deseq2.results.txt
    """
}
