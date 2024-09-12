from typing import List, Optional

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import (
    Aligner,
    BarcodeKit,
    Protocol,
    SampleSheet,
    initialize,
    nextflow_runtime,
)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_nanoseq(
    run_name: str,
    input: List[SampleSheet],
    email: Optional[str],
    multiqc_title: Optional[str],
    input_path_source: str,
    input_path_file: Optional[LatchFile],
    input_path_dir: Optional[LatchDir],
    barcode_kit: Optional[BarcodeKit],
    barcode_both_ends: bool,
    trim_barcodes: bool,
    gpu_cluster_options: Optional[str],
    qcat_detect_middle: bool,
    skip_demultiplexing: bool,
    run_nanolyse: bool,
    nanolyse_fasta: Optional[LatchFile],
    stranded: bool,
    save_align_intermeds: bool,
    skip_alignment: bool,
    call_variants: bool,
    split_mnps: bool,
    deepvariant_gpu: bool,
    phase_vcf: bool,
    skip_vc: bool,
    skip_sv: bool,
    skip_quantification: bool,
    skip_differential_analysis: bool,
    skip_fusion_analysis: bool,
    skip_modification_analysis: bool,
    skip_xpore: bool,
    skip_m6anet: bool,
    skip_bigbed: bool,
    skip_bigwig: bool,
    skip_nanoplot: bool,
    skip_fastqc: bool,
    skip_multiqc: bool,
    skip_qc: bool,
    protocol: Protocol = Protocol.cDNA,
    gpu_device: Optional[str] = "auto",
    qcat_min_score: Optional[int] = 60,
    aligner: Aligner = Aligner.minimap2,
    variant_caller: Optional[str] = "medaka",
    structural_variant_caller: Optional[str] = "sniffles",
    quantification_method: Optional[str] = "bambu",
    jaffal_ref_dir: Optional[str] = str("for_jaffal"),
    outdir: LatchOutputDir = LatchOutputDir("latch:///Nanoseq"),
) -> None:
    """
    nfcore/nanoseq is a bioinformatics analysis pipeline for Nanopore DNA/RNA sequencing data that can be used to perform basecalling, demultiplexing, QC, alignment, and downstream analysis.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3697959-1073c8)](https://doi.org/10.5281/zenodo.3697959)

    **nfcore/nanoseq** is a bioinformatics analysis pipeline for Nanopore DNA/RNA sequencing data that can be used to perform basecalling, demultiplexing, QC, alignment, and downstream analysis.

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    ## Pipeline Summary

    1. Demultiplexing ([`qcat`](https://github.com/nanoporetech/qcat); _optional_)
    2. Raw read cleaning ([NanoLyse](https://github.com/wdecoster/nanolyse); _optional_)
    3. Raw read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot), [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    4. Alignment ([`GraphMap2`](https://github.com/lbcb-sci/graphmap2) or [`minimap2`](https://github.com/lh3/minimap2))
       - Both aligners are capable of performing unspliced and spliced alignment. Sensible defaults will be applied automatically based on a combination of the input data and user-specified parameters
       - Each sample can be mapped to its own reference genome if multiplexed in this way
       - Convert SAM to co-ordinate sorted BAM and obtain mapping metrics ([`samtools`](http://www.htslib.org/doc/samtools.html))
    5. Create bigWig ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/)) and bigBed ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedToBigBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)) coverage tracks for visualisation
    6. DNA specific downstream analysis:
       - Short variant calling ([`medaka`](https://github.com/nanoporetech/medaka), [`deepvariant`](https://github.com/google/deepvariant), or [`pepper_margin_deepvariant`](https://github.com/kishwarshafin/pepper))
       - Structural variant calling ([`sniffles`](https://github.com/fritzsedlazeck/Sniffles) or [`cutesv`](https://github.com/tjiangHIT/cuteSV))
    7. RNA specific downstream analysis:
       - Transcript reconstruction and quantification ([`bambu`](https://bioconductor.org/packages/release/bioc/html/bambu.html) or [`StringTie2`](https://ccb.jhu.edu/software/stringtie/))
         - bambu performs both transcript reconstruction and quantification
         - When StringTie2 is chosen, each sample can be processed individually and combined. After which, [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/) will be used for both gene and transcript quantification.
       - Differential expression analysis ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and/or [`DEXSeq`](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html))
       - RNA modification detection ([`xpore`](https://github.com/GoekeLab/xpore) and/or [`m6anet`](https://github.com/GoekeLab/m6anet))
       - RNA fusion detection ([`JAFFAL`](https://github.com/Oshlack/JAFFA))
    8. Present QC for raw read and alignment results ([`MultiQC`](https://multiqc.info/docs/))

    ## Documentation

    The nf-core/nanoseq pipeline comes with documentation about the pipeline [usage](https://nf-co.re/nanoseq/usage), [parameters](https://nf-co.re/nanoseq/parameters) and [output](https://nf-co.re/nanoseq/output).

    ## Credits

    nf-core/nanoseq was originally written by [Chelsea Sawyer](https://github.com/csawye01) and [Harshil Patel](https://github.com/drpatelh) from [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) for use at [The Francis Crick Institute](https://www.crick.ac.uk/), London. Other primary contributors include [Laura Wratten](https://github.com/lwratten), [Ying Chen](https://github.com/cying111), [Yuk Kei Wan](https://github.com/yuukiiwa) and [Jonathan Goeke](https://github.com/jonathangoeke) from the [Genome Institute of Singapore](https://www.a-star.edu.sg/gis), [Christopher Hakkaart](https://github.com/christopher-hakkaart) from [Institute of Medical Genetics and Applied Genomics](https://www.medizin.uni-tuebingen.de/de/das-klinikum/einrichtungen/institute/medizinische-genetik-und-angewandte-genomik), Germany, and [Johannes Alneberg](https://github.com/alneberg) and [Franziska Bonath](https://github.com/FranBonath) from [SciLifeLab](https://www.scilifelab.se/), Sweden.

    Many thanks to others who have helped out along the way too, including (but not limited to): [@crickbabs](https://github.com/crickbabs), [@AnnaSyme](https://github.com/AnnaSyme), [@ekushele](https://github.com/ekushele).

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/nanoseq) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
        run_name=run_name,
        pvc_name=pvc_name,
        input=input,
        protocol=protocol,
        outdir=outdir,
        email=email,
        multiqc_title=multiqc_title,
        input_path_source=input_path_source,
        input_path_file=input_path_file,
        input_path_dir=input_path_dir,
        barcode_kit=barcode_kit,
        barcode_both_ends=barcode_both_ends,
        trim_barcodes=trim_barcodes,
        gpu_device=gpu_device,
        gpu_cluster_options=gpu_cluster_options,
        qcat_min_score=qcat_min_score,
        qcat_detect_middle=qcat_detect_middle,
        skip_demultiplexing=skip_demultiplexing,
        run_nanolyse=run_nanolyse,
        nanolyse_fasta=nanolyse_fasta,
        aligner=aligner,
        stranded=stranded,
        save_align_intermeds=save_align_intermeds,
        skip_alignment=skip_alignment,
        call_variants=call_variants,
        variant_caller=variant_caller,
        structural_variant_caller=structural_variant_caller,
        split_mnps=split_mnps,
        deepvariant_gpu=deepvariant_gpu,
        phase_vcf=phase_vcf,
        skip_vc=skip_vc,
        skip_sv=skip_sv,
        quantification_method=quantification_method,
        skip_quantification=skip_quantification,
        skip_differential_analysis=skip_differential_analysis,
        jaffal_ref_dir=jaffal_ref_dir,
        skip_fusion_analysis=skip_fusion_analysis,
        skip_modification_analysis=skip_modification_analysis,
        skip_xpore=skip_xpore,
        skip_m6anet=skip_m6anet,
        skip_bigbed=skip_bigbed,
        skip_bigwig=skip_bigwig,
        skip_nanoplot=skip_nanoplot,
        skip_fastqc=skip_fastqc,
        skip_multiqc=skip_multiqc,
        skip_qc=skip_qc,
    )


LaunchPlan(
    nf_nf_core_nanoseq,
    "Test Data",
    {
        "input": [
            SampleSheet(
                group="sample",
                replicate=1,
                barcode=6,
                input_file=None,
                fasta=LatchFile(
                    "s3://latch-public/nf-core/nanoseq/test_data/hg19_KCMF1.fa"
                ),
                gtf=None,
            ),
            SampleSheet(
                group="sample",
                replicate=2,
                barcode=12,
                input_file=None,
                fasta=LatchFile(
                    "s3://latch-public/nf-core/nanoseq/test_data/hg19_KCMF1.fa"
                ),
                gtf=None,
            ),
        ],
        "run_name": "Test_Run",
        "protocol": Protocol.DNA,
        "input_path_file": LatchFile(
            "s3://latch-public/nf-core/nanoseq/test_data/sample_nobc_dx.fastq.gz"
        ),
        "barcode_kit": BarcodeKit.NBD103_104,
        "skip_quantification": True,
        "skip_bigwig": True,
        "skip_bigbed": True,
        "skip_fusion_analysis": True,
        "skip_modification_analysis": True,
    },
)
