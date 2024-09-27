process SYLPH_SKETCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/sylph:0.6.1--h4ac6f70_0' :
        'quay.io/biocontainers/sylph@sha256:781d32f6e29a5ef8140f0f459f645525f43d021d0b5388e6caf2071d2e33ffd4' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path('*.sylsp'), optional: true, emit: sketch
    tuple val(meta), path('*.syldb'), optional: true, emit: genome_sketch
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq = input.getName().endsWith("fastq.gz") ? "-r $input" : ""
    def fasta = input.getName().endsWith("fasta") ? "-g $input" : ""

    """
    sylph sketch \\
      $args \\
      $fastq \\
      $fasta \\
      -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$( sylph -V 2>&1 | sed 's/sylph //' ))
    END_VERSIONS
    """
}
