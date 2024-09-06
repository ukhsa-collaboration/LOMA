process ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/biopython@sha256:f85da084589c559a5b99f3075cd4c5eeb1d846d37bc43dd0ee9d1b834dbafbbb' :
        'quay.io/biocontainers/biopython@sha256:f85da084589c559a5b99f3075cd4c5eeb1d846d37bc43dd0ee9d1b834dbafbbb' }"

    input:
    tuple val(meta), path(fasta, stageAs: "input_bins/*")

    output:
    tuple val(meta), path("*.assembly_stats.csv"), emit: asm_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assembly_stats.py \\
         --fasta_dir input_bins/ \\
         --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assembly_stats.py: \$(assembly_stats.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}
