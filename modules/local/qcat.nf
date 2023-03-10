process QCAT {
    tag "$input_path"
    label 'process_medium'

    conda "bioconda::qcat=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qcat:1.1.0--py_0' :
        'quay.io/biocontainers/qcat:1.1.0--py_0' }"

    input:
    path input_path

    output:
    path "fastq/*.fastq.gz", emit: fastq
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def detect_middle = params.qcat_detect_middle ? "--detect-middle $params.qcat_detect_middle" : ""
    """
    ## Unzip fastq file
    ## qcat doesnt support zipped files yet
    FILE=$input_path
    if [[ \$FILE == *.gz ]]
    then
    zcat $input_path > unzipped.fastq
    FILE=unzipped.fastq
    fi
    qcat  \\
    -f \$FILE \\
    -b ./fastq \\
    --kit $params.barcode_kit \\
    --min-score $params.qcat_min_score \\
    $detect_middle

    ## Zip fastq files (cannot find pigz command)
    gzip fastq/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qcat: \$(qcat --version 2>&1 | sed 's/^.*qcat //; s/ .*\$//')
    END_VERSIONS
    """
}
