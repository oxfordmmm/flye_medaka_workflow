process SUBSAMPLE {
    tag {sample}

    label 'short'

    input:
    tuple val(sample), path('reads.fastq.gz')

    output:
    tuple val(sample), path("${sample}.fastq.gz"), emit: fq

    script:
    """
    rasusa reads reads.fastq.gz --coverage 500 --genome-size 2.3m | gzip > ${sample}.fastq.gz
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
    label 'gpu'
    tag {sample}
    cpus 2


    label 'short'

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"}

    input:
    tuple val(sample), path('contigs.fasta'), path('reads.fastq.gz')
    //each path('model.tar.gz')

    output:
    tuple val(sample), path('output/consensus.fasta'), emit: fasta

    script:
    """
    medaka_consensus \
	-i reads.fastq.gz \
	-d contigs.fasta \
	-o output \
	-t ${task.cpus} 
    """
    stub:
    """
    mkdir output
    touch output/consensus.fasta
    """
}

process DNAAPLER {
    tag {sample}
    cpus 2

    publishDir "assemblies/${task.process.replaceAll(":","_")}", mode: 'copy', saveAs: { filename -> "${sample}.fasta"} 

    input:
    tuple val(sample),  path('input_assembly.fasta')
    
    output:
    tuple val(sample), path("dnaapler_out/${sample}*.fasta"), emit: fasta

    script:
    """
    dnaapler all \
        -i input_assembly.fasta -o dnaapler_out -p ${sample} -t ${task.cpus}


    """
    stub:
    """
    mkdir assembly
    touch assembly/assembly.fasta
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
    errorStrategy 'ignore'
    
    publishDir "bakta/", mode: 'copy'

    cpus 1

    input:
        tuple val(sample), path("consensus.fa"), path('db')

    output:
        tuple val(sample), path("${sample}_bakta/${sample}.gff3"), emit: gff3
        tuple val(sample), path("${sample}_bakta"), emit: fol


    script:
    """
    bakta --db db/db-light \
        -t ${task.cpus} \
        --prefix ${sample} \
        -o ${sample}_bakta/ \
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
    errorStrategy 'ignore'

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
            --isolate \
            -o ${sample}
    """
    stub:
    """
    mkdir ${sample}
    touch ${sample}/contigs.fasta
    """
}

process MLST {
    tag {sample}

    label 'short'

    publishDir "MLST/${task.process.replaceAll(":","_")}", mode: 'copy'

    input:
    tuple val(sample), path("${sample}.fasta")

    output:
    path("${sample}.tsv"), optional: true

    script:
    if (params.mlst_scheme != '') {
        """
        mlst --scheme $params.mlst_scheme --legacy ${sample}.fasta > ${sample}.tsv
        """
    } else {
        """
        mlst ${sample}.fasta > ${sample}.tsv
        """
    }
    stub:
    """
    touch ${sample}.tsv
    """
}

process COMBINE_MLST {
    label 'pandas'

    publishDir "MLST/${task.process.replaceAll(":","_")}", mode: 'copy'

    input:
    path('MLST???.tsv')
    each path('meta.csv')

    output:
    path("MLST_REPORT.csv")

    script:
    """
    combine_mlst.py -i MLST*.tsv -m meta.csv -o MLST_REPORT.csv -s $params.mlst_scheme
    """
    stub:
    """
    touch meta.csv
    """
}

process DNADIFF {
    tag {sample1 + ' ' + sample2}

    publishDir "DNAdiff/${task.process.replaceAll(":","_")}/${sample1}/${sample2}/", mode: 'copy'

    input:
    tuple val(sample1), path("${sample1}.fasta"), val(sample2), path("${sample2}.fasta")

    output:
    tuple val(sample1), val(sample2), path("out.report"), optional: true

    script:
    """
    dnadiff ${sample1}.fasta ${sample2}.fasta 
    """
    stub:
    """
    touch ${sample1}_${sample2}.report
    """
}

process PARSE_DD_REPORT {
    tag {sample1 + ' ' + sample2}
    label 'pandas'

    publishDir "DD_REPORTS/${task.process.replaceAll(":","_")}", mode: 'copy'

    input:
    tuple val(sample1), val(sample2), path("out.report")

    output:
    path("${sample1}_${sample2}.csv")

    script:
    """
    parse_dnadiff_report.py --sample_name1 ${sample1} \
        --sample_name2 ${sample2} \
        -i out.report \
        -o ${sample1}_${sample2}.csv 
    """
    stub:
    """
    touch ${sample1}_${sample2}.csv
    """
}

process SNP_MATRIX {
    label 'pandas'

    publishDir "SNP_MATRIX/${task.process.replaceAll(":","_")}", mode: 'copy'

    input:
    path("SNPs??????.csv")

    output:
    path("SNP_matrix.csv")

    script:
    """
    matrix.py -i SNPs*.csv -o SNP_matrix.csv
    """
    stub:
    """
    touch SNP_matrix.csv
    """
}

process FILTER_CONTIGS {
    tag {sample1}

    publishDir "Filtered_contigs/${task.process.replaceAll(":","_")}/", mode: 'copy'

    input:
    tuple val(sample1), path("${sample1}.fasta")

    output:
    tuple val(sample1), path("*${sample1}_filt.fasta"), optional: true

    script:
    """
    filter_contigs.py ${sample1}.fasta ${sample1}_filt.fasta 1500000 ${sample1}
    """
    stub:
    """
    touch ${sample1}_filt.fasta
    """
}