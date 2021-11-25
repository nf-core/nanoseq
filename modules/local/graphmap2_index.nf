// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GRAPHMAP2_INDEX {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"

    input:
    tuple path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), val(annotation_str)

    output:
    tuple path(fasta), path(sizes), path(gtf), val(bed), val(is_transcripts), path("*.gmidx"), val(annotation_str), emit: index
    path "versions.yml" , emit: versions

    script:
    def preset = (params.protocol == 'DNA' || is_transcripts) ? "" : "-x rnaseq"
    // TODO pipeline: Should be staging gtf file properly as an input
    def junctions = (params.protocol != 'DNA' && !is_transcripts && gtf) ? "--gtf $gtf" : ""
    """
    graphmap2 \\
        align \\
        $preset \\
        $junctions \\
        -t $task.cpus \\
        -I \\
        -r $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
