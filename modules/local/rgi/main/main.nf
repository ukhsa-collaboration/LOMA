process RGI_MAIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_1' :
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    tuple val(meta), path("temp/")  , emit: tmp
    env VER                        , emit: tool_version
    env DBVER                      , emit: db_version
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mv $fasta ${meta.id}

    rgi \\
        main \\
        $args \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence ${meta.id}

    mkdir temp/
    mv *.xml *.fsa *.{nhr,nin,nsq} *.draft *.potentialGenes *{variant,rrna,protein,predictedGenes,overexpression,homolog}.json temp/

    VER=\$(rgi main --version)
    DBVER=\$(rgi database --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(echo \$VER)
        rgi-database: \$(echo \$DBVER)
    END_VERSIONS
    """
}
