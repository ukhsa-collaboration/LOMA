process FILTER_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/samtools:1.20--h50ea8bc_0' :
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(sub_contigs), path(bam)

    output:
    tuple val(meta), path(sub_contigs), path("*.fastq.gz"), emit: sp_sub_fastq
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools faidx $sub_contigs
    awk '{print \$1,0,\$2+1}' ${sub_contigs}.fai > ${prefix}.sub_ref.bed
    samtools view -@ ${task.cpus} -h -L ${prefix}.sub_ref.bed $bam | samtools sort -n -@ 1 -o ${prefix}.sub_filt.bam -

    samtools fastq -@ ${task.cpus} -0 ${prefix}.${meta.bin}.fastq.gz ${prefix}.sub_filt.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
