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
include {ROARY} from './modules/assemble.nf'
include {SPADES} from './modules/assemble.nf'



workflow {
    // channels
    Channel.fromPath( "${params.inputFastq}/*" )
           .map{ file -> tuple(file.simpleName, file) }
           .set{ ch_fqs }

    Channel.fromPath("${params.inputmodel}")
        .set{ ch_model }

    Channel.fromFilePairs("illumina_fqs/*{1,2}.fastq.gz")
        .map{ row -> tuple(row[0], row[1][0], row[1][1])}
        .set{illumina_fqs}


    main:
    SUBSAMPLE(ch_fqs)

    ASSEMBLE(SUBSAMPLE.out.fq)

    SEQKIT(ASSEMBLE.out.fasta)

    POLISH(ASSEMBLE.out.fasta.combine(SUBSAMPLE.out.fq, by:0), ch_model)

    //PROKKA(ASSEMBLE.out.fasta)

    BAKTA_DOWNLOAD()

    SPADES(illumina_fqs)

    contigs=POLISH.out.fasta.mix(SPADES.out.fasta)

    BAKTA(contigs.combine(BAKTA_DOWNLOAD.out))

    gff3s=BAKTA.out.gff3
        .map{ row -> row[1]}
        .collect()

    ROARY(gff3s)

}