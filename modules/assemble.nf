process SUBSAMPLE {
    tag {sample}

    label 'short'

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

    label 'short'

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

    label 'short'

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), path('contigs.fasta'), path('reads.fastq.gz')
    each path('model.tar.gz')

    output:
    tuple val(sample), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus \
	-i reads.fastq.gz \
	-d contigs.fasta \
	-o output \
	-t ${task.cpus} \
	-m model.tar.gz
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}

process PROKKA {
    tag { sample }
    
    publishDir "prokka/", mode: 'copy'

    cpus 8

    input:
        tuple val(sample), path("consensus.fa")

    output:
        tuple val(sample), path("${sample}_prokka"), emit: ch_out

    script:
    """
    prokka --cpus ${task.cpus} \
        --addgenes \
        --force \
        --compliant \
        --outdir ${sample}_prokka \
        --prefix ${sample} consensus.fa
    """
    stub:
    """
    mkdir ${sample}_prokka
    touch ${sample}_prokka/${sample}.gbk 
    """
}

process BAKTA_DOWNLOAD {
    label 'online'

    output:
    path('db')

    script:
    """
    bakta_db download --output db --type light
    """

}

process BAKTA {
    tag { sample }
    
    publishDir "bakta/", mode: 'copy'

    cpus 8

    input:
        tuple val(sample), path("consensus.fa"), path('db')

    output:
        tuple val(sample), path("${sample}.gff3"), emit: gff3

    script:
    """
    bakta --db db/db-light \
        -t ${task.cpus} \
        --prefix ${sample} \
        consensus.fa
    """
    stub:
    """
    mkdir ${sample}_prokka
    touch ${sample}_prokka/${sample}.gbk 
    """
}

process SEQKIT {
    tag {sample}

    publishDir "CSVs/${task.process.replaceAll(":","_")}", mode: 'copy'

    input:
    tuple val(sample),  path('contigs.fasta')

    output:
    tuple val(sample),  path("${sample}.tsv")

    script:
    """
    seqkit stats -T --all contigs.fasta > ${sample}.tsv
    """
    stub:
    """
    touch ${sample}.tsv
    """
}

process ROARY {

    publishDir "roary/${task.process.replaceAll(":","_")}", mode: 'copy'

    cpus 8

    input:
    path('*')

    output:
    path("roary_output")

    script:
    """
    roary -p ${task.cpus} -f roary_output *.gff3
    """
    stub:
    """
    mkdir roary_output
    """
}

process SPADES {
    tag {sample}
    cpus 4

    label 'short'

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 

    input:
    //tuple val(sample), val(refname), val(reads)
    tuple val(sample),  file('reads1.fastq.gz'), file('reads2.fastq.gz')

    output:
    tuple val(sample), file("${sample}/contigs.fasta"), emit: fasta

    script:
    """
    spades.py -t ${task.cpus} \
            --pe1-1 reads1.fastq.gz \
            --pe1-2 reads2.fastq.gz \
            -o ${sample}
    """
    stub:
    """
    mkdir ${sample}
    touch ${sample}/contigs.fasta
    """
}