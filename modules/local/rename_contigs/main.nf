process RENAME_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(bins), path(unbinned)

    output:
    tuple val(meta), path("assembled_bins/*[!unbinned].fasta"), optional: true, emit: assembled_bins
    tuple val(meta), path("*.complete_assembly.fasta"), emit: complete_assembly
    tuple val(meta), path("assembled_bins/*.unbinned.fasta"), optional: true, emit: unbinned
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir assembled_bins
    rm -f empty.fa
    mv *.fa assembled_bins/

    rename_bins.py \\
       --input_dir assembled_bins/ \\
       --prefix ${meta.id}

    cat assembled_bins/* > ${prefix}.complete_assembly.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rename_bins.py: \$(rename_bins.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}
