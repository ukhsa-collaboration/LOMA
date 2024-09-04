process ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_single'

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
