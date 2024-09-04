process FETCH_SEQSUM {
    tag "$meta"
    label 'process_low'

    input:
    val(meta)

    output:
    tuple val(meta), path("sequencing_summary*"), optional: true, emit: sequencing_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    process_seqsum.py \\
       --input-dir $params.run_dir \\
       --run_id ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        process_seqsum.py: \$(process_seqsum.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}
