#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// include modules
include {SUBSAMPLE} from './modules/assemble.nf'
include {ASSEMBLE} from './modules/assemble.nf'
include {POLISH} from './modules/assemble.nf'


workflow {
    // channels
    Channel.fromPath( "${params.inputFastq}/*" )
           .map{ file -> tuple(file.simpleName, file) }
           .set{ ch_fqs }


    main:
    SUBSAMPLE(ch_fqs)

    ASSEMBLE(SUBSAMPLE.out.fq)

    POLISH(ASSEMBLE.out.fasta.combine(SUBSAMPLE.out.fq, by:0))

}