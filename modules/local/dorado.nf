process DORADO {
    label 'process_medium'

    container "docker.io/ontresearch/dorado"

    input:
    path(input_path)
    val meta
    path dorado_config
    path dorado_model

    output:
    path "*.fastq.gz"                    , emit: fastq
    path "versions.yml"                  , emit: versions

    script:
    def fast5_dir_path = workflow.profile.contains('test') ? "input_path" : "$input_path"
    def trim_barcodes = params.trim_barcodes ? "--trim_barcodes" : ""
    def barcode_kit  = params.barcode_kit ? "--barcode_kits $params.barcode_kit" : ""
    def barcode_ends = params.barcode_both_ends ? "--require_barcodes_both_ends" : ""
    def proc_options = params.dorado_gpu ? "--device $params.gpu_device --num_callers $task.cpus --cpu_threads_per_caller $params.dorado_cpu_threads --gpu_runners_per_device $params.dorado_gpu_runners" : "--num_callers 2 --cpu_threads_per_caller ${task.cpus/2}"
    def config   = "--flowcell $params.flowcell --kit $params.kit"
    if (params.dorado_config) config = file(params.dorado_config).exists() ? "--config ./$dorado_config" : "--config $params.dorado_config"
    def model    = ""
    if (params.dorado_model)  model  = file(params.dorado_model).exists() ? "--model ./$dorado_model" : "--model $params.dorado_model"
    """
    dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.1.0
    dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.1.0 $input_path --device cpu --emit-fastq > basecall.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(echo \$(dorado --version 2>&1) | sed -r 's/.{81}//')
    END_VERSIONS

    gzip basecall.fastq
    """
}

