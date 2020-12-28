// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GRAPHMAP2_ALIGN {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"

    input:
    tuple val(sample), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)
    
    output:
    tuple val(sample), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
    // TODO pipeline: Should be staging gtf file properly as an input
    def junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf $gtf" : ""
    """
    graphmap2 \\
        align \\
        $preset \\
        $junctions \\
        -t $task.cpus \\
        -r $fasta \\
        -i $index \\
        -d $fastq \\
        -o ${sample}.sam \\
        --extcigar
    """
}
