workflow DORADO {
    take:
    path samplesheet
    path input_files //path to dir containing fastqs/ubams
    path reference_fasta
    // path reference_fai
    val barcode_kit

    main:
    //TODO make demux optional
    DORADO_DEMUX(samplesheet, fastq, barcode_kit) //must output fastqs!
    sample_fastq = DORADO_DEMUX.out.fastq //{[meta, files]}

    CHOPPER(sample_fastq)
    sample_clean_fastq = CHOPPER.out.clean_fastq

    DORADO_ALIGN(sample_fastq, reference_fasta)
    sample_bams = DORADO_ALIGN.out.bam //{[meta1, bam1], [meta2, bam2]}

    // SAMTOOLS_SORT(sample_bams)
    // sample_sorted_bams = SAMTOOLS_SORT.out.bam

    // SAMTOOLS_INDEX(sample_sorted_bams)
    // sample_bai = SAMTOOLS_INDEX.out.bai //{[meta1, bai1], [meta2, bai2]}

    // sample_output_bams = sample_sorted_bams.join(sample_bai) //{[meta1, bam1, bai1], [meta2, bam2, bai2]}

    emit:
    sample_bams

}
