process SUBSAMPLE {
    tag {sample}

    input:
    tuple val(sample), path('reads.fastq.gz')

    output:
    tuple val(sample), path("${sample}.fastq.gz"), emit: fq

    script:
    """
    rasusa --input reads.fastq.gz --coverage 100 --genome-size 2.4m | gzip > ${sample}.fastq.gz
    """ 
    stub:
    """
    touch ${sample}.fastq.gz
    """
}

process ASSEMBLE {
    tag {sample}
    cpus 4

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 

    input:
    //tuple val(sample), val(refname), val(reads)
    tuple val(sample),  file('reads.fastq.gz')

    output:
    tuple val(sample), file('assembly/assembly.fasta'), emit: fasta

    script:
    """
    flye -o assembly \
	--plasmids \
	--meta \
	--threads ${task.cpus} \
	--nano-hq reads.fastq.gz
    """
    stub:
    """
    mkdir assembly
    touch assembly/assembly.fasta
    """
}

process POLISH {
    label 'medaka'
    tag {sample}
    cpus 2

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), path('contigs.fasta'), path('reads.fastq.gz'), val(model)

    output:
    tuple val(sample), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus \
	-i reads.fastq.gz \
	-d contigs.fasta \
	-o output \
	-t ${task.cpus} \
	-m ${model}
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}