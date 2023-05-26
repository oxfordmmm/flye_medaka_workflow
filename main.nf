#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// include modules
include {SUBSAMPLE} from './modules/assemble.nf'
include {ASSEMBLE} from './modules/assemble.nf'
include {POLISH} from './modules/assemble.nf'
include {SEQKIT} from './modules/assemble.nf'
include {PROKKA} from './modules/assemble.nf'
include {BAKTA_DOWNLOAD} from './modules/assemble.nf'
include {BAKTA} from './modules/assemble.nf'


workflow {
    // channels
    Channel.fromPath( "${params.inputFastq}/*" )
           .map{ file -> tuple(file.simpleName, file) }
           .set{ ch_fqs }

    Channel.fromPath("${params.inputmodel}")
        .set{ ch_model }


    main:
    SUBSAMPLE(ch_fqs)

    ASSEMBLE(SUBSAMPLE.out.fq)

    SEQKIT(ASSEMBLE.out.fasta)

    POLISH(ASSEMBLE.out.fasta.combine(SUBSAMPLE.out.fq, by:0))

    //PROKKA(ASSEMBLE.out.fasta)

    BAKTA_DOWNLOAD()

    BAKTA(ASSEMBLE.out.fasta.combine(BAKTA_DOWNLOAD.out))


}