// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options       = [:]
def options          = initOptions(params.options)

process GET_NANOLYSE_FASTA {
    output:    
    path "*fasta.gz"  , emit: ch_nanolyse_fasta

    script:
    """
    wget https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz
    """
}
