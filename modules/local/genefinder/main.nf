process GENEFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/genefdep:latest' :
        'quay.io/djberger/genefdep:latest' }"

    input:
    tuple val(meta), path(assembly), path(reads)
    path(dbdir)

    output:
    tuple val(meta), path("*.xml"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def DB_name = meta.gene_DB

    """
    mkdir -p rundir
    mkdir -p outdir

    ONT_gene_finder.py \\
    $args \\
    -r rundir/ \\
    --file $reads \\
    -o outdir/ \\
    -db $dbdir/$DB_name/reference.fasta \\
    --threads 4 \\
    -x

    mv outdir/* ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ONT_gene_finder.py: \$(grep "^__version__" /data/PROJECTS/DB_ONT_PIPELINE_TEST/loma_v0.2/bin/ONT_gene_finder.py | cut -f3 -d " " | sed 's/"//g')
    END_VERSIONS
    """
}
