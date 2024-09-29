process RACON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/racon:1.4.20--h9a82719_1' :
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(meta), path(reads), path(assembly), path(paf)
    val f_round

    output:
    tuple val(meta), path(reads), path('*_assembly_consensus.fasta') , emit: improved_assembly
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def last_round = f_round ? true : false

    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        "${paf}" \\
        $args \\
        "${assembly}" > \\
        ${prefix}_assembly_consensus.fasta

    if [ "$last_round" == "true" ]; then
        mv ${prefix}_assembly_consensus.fasta ${prefix}_assembly_consensus.fasta.1
        cut -f1 -d " " ${prefix}_assembly_consensus.fasta.1 > ${prefix}_assembly_consensus.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        racon: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
