process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0' :
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(meta) , path(fasta)
    path(dbdir)

    output:
    tuple val(meta), path('workflow.txt'), path('*blastn*'), emit: txt, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def DB_name = meta.gene_DB

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`

    cp $dbdir/$DB_name/workflow.txt ./

    blastn \\
        -num_threads ${task.cpus} \\
        -db $dbdir/$DB_name/reference.fasta \\
        -query ${fasta_name} \\
        ${args} \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
