#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// define the input parameters
params.mlst_scheme = ''
params.external_fasta = ''

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
include {MLST} from './modules/assemble.nf'
include {COMBINE_MLST} from './modules/assemble.nf'
include {DNADIFF} from './modules/assemble.nf'
include {PARSE_DD_REPORT} from './modules/assemble.nf'
include {SNP_MATRIX} from './modules/assemble.nf'
include {FILTER_CONTIGS} from './modules/assemble.nf'
include {DNAAPLER} from './modules/assemble.nf'



workflow {
    // channels
    Channel.fromPath( "${params.inputFastq}/*" )
           .map{ file -> tuple(file.simpleName, file) }
           .set{ ch_fqs }

    Channel.fromPath("${params.inputmodel}")
        .set{ ch_model }

    Channel.fromPath("${params.meta}")
        .set{ ch_meta }

    Channel.fromFilePairs("illumina_fqs/*{1,2}.fastq.gz")
        .map{ row -> tuple(row[0], row[1][0], row[1][1])}
        .set{illumina_fqs}

    Channel.fromPath("${params.external_fasta}/*")
        .map{ file -> tuple(file.simpleName, file) }
        .set{ ch_external_fasta }


    main:
    SUBSAMPLE(ch_fqs)

    ASSEMBLE(SUBSAMPLE.out.fq)

    SEQKIT(ASSEMBLE.out.fasta)

    //POLISH(ASSEMBLE.out.fasta.combine(SUBSAMPLE.out.fq, by:0), ch_model)
    POLISH(ASSEMBLE.out.fasta.combine(SUBSAMPLE.out.fq, by:0))

    //PROKKA(ASSEMBLE.out.fasta)

    BAKTA_DOWNLOAD()

    SPADES(illumina_fqs)

    polish_contigs=POLISH.out.fasta.mix(SPADES.out.fasta)

    DNAAPLER(polish_contigs)

    contigs=DNAAPLER.out.mix(ch_external_fasta)

    FILTER_CONTIGS(contigs)

    //BAKTA(FILTER_CONTIGS.out.combine(BAKTA_DOWNLOAD.out))

    //gff3s=BAKTA.out.gff3
    //    .map{ row -> row[1]}
    //    .collect()

    //ROARY(gff3s)

    MLST(FILTER_CONTIGS.out)
    COMBINE_MLST(MLST.out.collect(), ch_meta)

    // paiwise comparison of each sample
    pairs=FILTER_CONTIGS.out.combine(FILTER_CONTIGS.out)
        .filter{ it[0] != it[2] }

    DNADIFF(pairs)

    // parse the report
    PARSE_DD_REPORT(DNADIFF.out)
    SNP_MATRIX(PARSE_DD_REPORT.out.collect())
}