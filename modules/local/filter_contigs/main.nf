process FILTER_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/samtools:1.20--h50ea8bc_0' :
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(sub_contigs), path(reference), path(bam), path(bai)

    output:
    tuple val(meta), path("*.sub_ref.fa"), path("*.sub_filt.bam"), path("*.bai"), optional: true, emit: filt_ref_bam
    path  "versions.yml"           , optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = reference.getName().endsWith(".gz") ? true : false
    def reference_name = reference.getName().replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $reference > $reference_name
    fi

    samtools faidx -r $sub_contigs $reference_name > ${prefix}.sub_ref.fa
    samtools faidx ${prefix}.sub_ref.fa
    awk '{print \$1,0,\$2+1}' ${prefix}.sub_ref.fa.fai > ${prefix}.sub_ref.bed
    samtools view -@ ${task.cpus} -L ${prefix}.sub_ref.bed -o ${prefix}.sub_filt.bam $bam

    samtools \\
        index \\
        -@ ${task.cpus} \\
        $args \\
        ${prefix}.sub_filt.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
