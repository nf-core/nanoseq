from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_gib": 100,
        }
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: LatchFile, protocol: str, outdir: typing.Optional[typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})]], email: typing.Optional[str], multiqc_title: typing.Optional[str], input_path: typing.Optional[LatchFile], barcode_kit: typing.Optional[str], barcode_both_ends: typing.Optional[bool], trim_barcodes: typing.Optional[bool], gpu_cluster_options: typing.Optional[str], qcat_detect_middle: typing.Optional[bool], skip_demultiplexing: typing.Optional[bool], run_nanolyse: typing.Optional[bool], nanolyse_fasta: typing.Optional[str], stranded: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], skip_alignment: typing.Optional[bool], call_variants: typing.Optional[bool], split_mnps: typing.Optional[bool], deepvariant_gpu: typing.Optional[bool], phase_vcf: typing.Optional[bool], skip_vc: typing.Optional[bool], skip_sv: typing.Optional[bool], skip_quantification: typing.Optional[bool], skip_differential_analysis: typing.Optional[bool], skip_fusion_analysis: typing.Optional[bool], skip_modification_analysis: typing.Optional[bool], skip_xpore: typing.Optional[bool], skip_m6anet: typing.Optional[bool], skip_bigbed: typing.Optional[bool], skip_bigwig: typing.Optional[bool], skip_nanoplot: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], gpu_device: typing.Optional[str], qcat_min_score: typing.Optional[int], aligner: typing.Optional[str], variant_caller: typing.Optional[str], structural_variant_caller: typing.Optional[str], quantification_method: typing.Optional[str], jaffal_ref_dir: typing.Optional[str]) -> None:
    try:
        shared_dir = Path("/nf-workdir")



        ignore_list = [
            "latch",
            ".latch",
            "nextflow",
            ".nextflow",
            "work",
            "results",
            "miniconda",
            "anaconda3",
            "mambaforge",
        ]

        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
            "-c",
            "latch.config",
                *get_flag('input', input),
                *get_flag('protocol', protocol),
                *get_flag('outdir', outdir),
                *get_flag('email', email),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('input_path', input_path),
                *get_flag('barcode_kit', barcode_kit),
                *get_flag('barcode_both_ends', barcode_both_ends),
                *get_flag('trim_barcodes', trim_barcodes),
                *get_flag('gpu_device', gpu_device),
                *get_flag('gpu_cluster_options', gpu_cluster_options),
                *get_flag('qcat_min_score', qcat_min_score),
                *get_flag('qcat_detect_middle', qcat_detect_middle),
                *get_flag('skip_demultiplexing', skip_demultiplexing),
                *get_flag('run_nanolyse', run_nanolyse),
                *get_flag('nanolyse_fasta', nanolyse_fasta),
                *get_flag('aligner', aligner),
                *get_flag('stranded', stranded),
                *get_flag('save_align_intermeds', save_align_intermeds),
                *get_flag('skip_alignment', skip_alignment),
                *get_flag('call_variants', call_variants),
                *get_flag('variant_caller', variant_caller),
                *get_flag('structural_variant_caller', structural_variant_caller),
                *get_flag('split_mnps', split_mnps),
                *get_flag('deepvariant_gpu', deepvariant_gpu),
                *get_flag('phase_vcf', phase_vcf),
                *get_flag('skip_vc', skip_vc),
                *get_flag('skip_sv', skip_sv),
                *get_flag('quantification_method', quantification_method),
                *get_flag('skip_quantification', skip_quantification),
                *get_flag('skip_differential_analysis', skip_differential_analysis),
                *get_flag('jaffal_ref_dir', jaffal_ref_dir),
                *get_flag('skip_fusion_analysis', skip_fusion_analysis),
                *get_flag('skip_modification_analysis', skip_modification_analysis),
                *get_flag('skip_xpore', skip_xpore),
                *get_flag('skip_m6anet', skip_m6anet),
                *get_flag('skip_bigbed', skip_bigbed),
                *get_flag('skip_bigwig', skip_bigwig),
                *get_flag('skip_nanoplot', skip_nanoplot),
                *get_flag('skip_fastqc', skip_fastqc),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('skip_qc', skip_qc)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
        print(flush=True)

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_nanoseq", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def nf_nf_core_nanoseq(input: LatchFile, protocol: str, outdir: typing.Optional[typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})]], email: typing.Optional[str], multiqc_title: typing.Optional[str], input_path: typing.Optional[LatchFile], barcode_kit: typing.Optional[str], barcode_both_ends: typing.Optional[bool], trim_barcodes: typing.Optional[bool], gpu_cluster_options: typing.Optional[str], qcat_detect_middle: typing.Optional[bool], skip_demultiplexing: typing.Optional[bool], run_nanolyse: typing.Optional[bool], nanolyse_fasta: typing.Optional[str], stranded: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], skip_alignment: typing.Optional[bool], call_variants: typing.Optional[bool], split_mnps: typing.Optional[bool], deepvariant_gpu: typing.Optional[bool], phase_vcf: typing.Optional[bool], skip_vc: typing.Optional[bool], skip_sv: typing.Optional[bool], skip_quantification: typing.Optional[bool], skip_differential_analysis: typing.Optional[bool], skip_fusion_analysis: typing.Optional[bool], skip_modification_analysis: typing.Optional[bool], skip_xpore: typing.Optional[bool], skip_m6anet: typing.Optional[bool], skip_bigbed: typing.Optional[bool], skip_bigwig: typing.Optional[bool], skip_nanoplot: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], gpu_device: typing.Optional[str] = 'auto', qcat_min_score: typing.Optional[int] = 60, aligner: typing.Optional[str] = 'minimap2', variant_caller: typing.Optional[str] = 'medaka', structural_variant_caller: typing.Optional[str] = 'sniffles', quantification_method: typing.Optional[str] = 'bambu', jaffal_ref_dir: typing.Optional[str] = 'for_jaffal') -> None:
    """
    nf-core/nanoseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, protocol=protocol, outdir=outdir, email=email, multiqc_title=multiqc_title, input_path=input_path, barcode_kit=barcode_kit, barcode_both_ends=barcode_both_ends, trim_barcodes=trim_barcodes, gpu_device=gpu_device, gpu_cluster_options=gpu_cluster_options, qcat_min_score=qcat_min_score, qcat_detect_middle=qcat_detect_middle, skip_demultiplexing=skip_demultiplexing, run_nanolyse=run_nanolyse, nanolyse_fasta=nanolyse_fasta, aligner=aligner, stranded=stranded, save_align_intermeds=save_align_intermeds, skip_alignment=skip_alignment, call_variants=call_variants, variant_caller=variant_caller, structural_variant_caller=structural_variant_caller, split_mnps=split_mnps, deepvariant_gpu=deepvariant_gpu, phase_vcf=phase_vcf, skip_vc=skip_vc, skip_sv=skip_sv, quantification_method=quantification_method, skip_quantification=skip_quantification, skip_differential_analysis=skip_differential_analysis, jaffal_ref_dir=jaffal_ref_dir, skip_fusion_analysis=skip_fusion_analysis, skip_modification_analysis=skip_modification_analysis, skip_xpore=skip_xpore, skip_m6anet=skip_m6anet, skip_bigbed=skip_bigbed, skip_bigwig=skip_bigwig, skip_nanoplot=skip_nanoplot, skip_fastqc=skip_fastqc, skip_multiqc=skip_multiqc, skip_qc=skip_qc)

