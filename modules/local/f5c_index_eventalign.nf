process F5C_INDEX_EVENTALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nanopolish==0.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/f5c:1.5--h56e2c18_1' :
        'quay.io/biocontainers/f5c:1.5--h56e2c18_1' }"

    input:
    tuple val(meta), path(genome), path(gtf), path(fastqgz), path(bam), path(bai)

    output:
    tuple val(meta), path(genome), path(gtf), path("*eventalign.txt"), path("*summary.txt"), emit: f5c_outputs
    path "versions.yml"                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    sample_summary = "$meta.id" +"_summary.txt"
    sample_eventalign = "$meta.id" +"_eventalign.txt"
    fastq="$meta.id"+".fastq"
    fast5 = "$meta.fast5"
    fastqi="$fastq"+"*"
    """
    gunzip -c $fastqgz > $fastq
    f5c index -d $fast5 $fastq
    echo $fastqi
    f5c eventalign --reads $fastq --bam $bam --genome $genome --scale-events --signal-index --rna --min-mapq 0 --min-recalib-events 1 --summary $sample_summary --threads $task.cpus > $sample_eventalign

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$( f5c -V | tail -c 4)
    END_VERSIONS
    """
}
