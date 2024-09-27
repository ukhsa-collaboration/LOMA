process SEQTK_FQCHK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/seqtk:1.4--he4a0461_1' : 
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(fastx)
    val(endseq_len)

    output:
    tuple val(meta), path("*.fw.txt"), path("*.rv.txt")     , emit: pbq_se
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqtk \\
        fqchk \\
        $args \\
        $fastx | awk '\$1<($endseq_len+1)' > ${prefix}.fw.txt

    seqtk \\
        seq \\
        -r \\
        $args \\
        $fastx | seqtk fqchk - | awk '\$1<($endseq_len+1)' > ${prefix}.rv.txt

    echo "done"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
