process KROCUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/krocus:1.0.3--pyhdfd78af_0' :
        'quay.io/biocontainers/krocus:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(fastq)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    organism_param_3 = meta.containsKey("krocus_scheme") ? "${meta.krocus_scheme}" : "Other"

    """
    if [ -e mlst_files ]; then
       rm -rf mlst_files
    fi

    krocus_database_downloader --species "${organism_param_3}"

    krocus --output ${prefix}.tsv mlst_files/ $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(krocus --version 2>&1))
    END_VERSIONS
    """
}
