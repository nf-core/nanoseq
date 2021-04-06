// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'meta.id') }

    conda     (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    container "quay.io/biocontainers/minimap2:2.17--hed695b0_3"

    input:
    tuple val(meta), path(fastq), path(fasta), path(sizes), val(gtf), val(bed), val(is_transcripts), path(index)
    
    output:
    tuple val(meta), path(sizes), val(is_transcripts), path("*.sam"), emit: align_sam

    script:
    def preset    = (params.protocol == 'DNA' || is_transcripts) ? "-ax map-ont" : "-ax splice"
    def kmer      = (params.protocol == 'directRNA') ? "-k14" : ""
    def stranded  = (params.stranded || params.protocol == 'directRNA') ? "-uf" : ""
    // TODO pipeline: Should be staging bed file properly as an input
    def junctions = (params.protocol != 'DNA' && bed) ? "--junc-bed ${file(bed)}" : ""
    """
    minimap2 \\
        $preset \\
        $kmer \\
        $stranded \\
        $junctions \\
        -t $task.cpus \\
        $index \\
        $fastq > ${meta.id}.sam
    """
}
