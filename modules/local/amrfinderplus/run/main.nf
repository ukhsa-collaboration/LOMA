process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ncbi/amr:3.12.8-2024-07-22.1':
        'ncbi/amr:3.12.8-2024-07-22.1' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
    def is_compressed_db = db.getName().endsWith(".gz") ? true : false
    prefix = task.ext.prefix ?: "${meta.id}"

    organism_param = !(meta.amfr.matches("")) ? "--organism ${meta.amfr} --mutation_all ${prefix}-mutations.tsv" : ""
    fasta_name = fasta.getName().replace(".gz", "")
    fasta_param = "-n"
    if (meta.containsKey("is_proteins")) {
        if (meta.is_proteins) {
            fasta_param = "-p"
        }
    }
    """
    amrfinder \\
        $fasta_param $fasta_name \\
        $organism_param \\
        $args \\
        --database $db/ \\
        --threads $task.cpus > ${prefix}.tsv


    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    amrfinderplus: \$(amrfinder --version)
    amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database $db --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    """
}
