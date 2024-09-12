import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Union

import requests
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)

sys.stdout.reconfigure(line_buffering=True)


@dataclass
class SampleSheet:
    group: str
    replicate: int
    barcode: Optional[int]
    input_file: Optional[LatchFile]
    fasta: Optional[Union[LatchFile, str]]
    gtf: Optional[LatchFile]


class Protocol(Enum):
    cDNA = "cDNA"
    DNA = "DNA"
    directRNA = "directRNA"


class BarcodeKit(Enum):
    Auto = "Auto"
    RBK001 = "RBK001"
    RBK004 = "RBK004"
    NBD103_104 = "NBD103/NBD104"
    NBD114 = "NBD114"
    NBD104_114 = "NBD104/NBD114"
    PBC001 = "PBC001"
    PBC096 = "PBC096"
    RPB004_RLB001 = "RPB004/RLB001"
    RPB004_LWB001 = "RPB004/LWB001"
    RAB204 = "RAB204"
    VMK001 = "VMK001"


class Aligner(Enum):
    minimap2 = "minimap2"
    graphmap2 = "graphmap2"


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        # "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


input_construct_samplesheet = metadata._nextflow_metadata.parameters[
    "input"
].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    run_name: str,
    pvc_name: str,
    input: List[SampleSheet],
    protocol: Protocol,
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
    aligner: Aligner,
    stranded: bool,
    save_align_intermeds: bool,
    skip_alignment: bool,
    call_variants: bool,
    variant_caller: Optional[str],
    structural_variant_caller: Optional[str],
    split_mnps: bool,
    deepvariant_gpu: bool,
    phase_vcf: bool,
    skip_vc: bool,
    skip_sv: bool,
    quantification_method: Optional[str],
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
    gpu_device: Optional[str],
    qcat_min_score: Optional[int],
    jaffal_ref_dir: Optional[str],
    outdir: LatchOutputDir,
) -> None:
    shared_dir = Path("/nf-workdir")

    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
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
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag("protocol", protocol),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("email", email),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("input_path", input_path_dir if input_path_dir else input_path_file),
        *get_flag("barcode_kit", barcode_kit),
        *get_flag("barcode_both_ends", barcode_both_ends),
        *get_flag("trim_barcodes", trim_barcodes),
        *get_flag("gpu_device", gpu_device),
        *get_flag("gpu_cluster_options", gpu_cluster_options),
        *get_flag("qcat_min_score", qcat_min_score),
        *get_flag("qcat_detect_middle", qcat_detect_middle),
        *get_flag("skip_demultiplexing", skip_demultiplexing),
        *get_flag("run_nanolyse", run_nanolyse),
        *get_flag("nanolyse_fasta", nanolyse_fasta),
        *get_flag("aligner", aligner),
        *get_flag("stranded", stranded),
        *get_flag("save_align_intermeds", save_align_intermeds),
        *get_flag("skip_alignment", skip_alignment),
        *get_flag("call_variants", call_variants),
        *get_flag("variant_caller", variant_caller),
        *get_flag("structural_variant_caller", structural_variant_caller),
        *get_flag("split_mnps", split_mnps),
        *get_flag("deepvariant_gpu", deepvariant_gpu),
        *get_flag("phase_vcf", phase_vcf),
        *get_flag("skip_vc", skip_vc),
        *get_flag("skip_sv", skip_sv),
        *get_flag("quantification_method", quantification_method),
        *get_flag("skip_quantification", skip_quantification),
        *get_flag("skip_differential_analysis", skip_differential_analysis),
        *get_flag("jaffal_ref_dir", jaffal_ref_dir),
        *get_flag("skip_fusion_analysis", skip_fusion_analysis),
        *get_flag("skip_modification_analysis", skip_modification_analysis),
        *get_flag("skip_xpore", skip_xpore),
        *get_flag("skip_m6anet", skip_m6anet),
        *get_flag("skip_bigbed", skip_bigbed),
        *get_flag("skip_bigwig", skip_bigwig),
        *get_flag("skip_nanoplot", skip_nanoplot),
        *get_flag("skip_fastqc", skip_fastqc),
        *get_flag("skip_multiqc", skip_multiqc),
        *get_flag("skip_qc", skip_qc),
    ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///your_log_dir/nf_nf_core_nanoseq", name, "nextflow.log"
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
