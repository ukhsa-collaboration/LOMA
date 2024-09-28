process SEMIBIN_SINGLEEASYBIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0' :
        'https://depot.galaxyproject.org/singularity/semibin:2.1.0--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("*.csv"), optional: true, emit: csv
    tuple val(meta), path("output_bins/*"), optional: true, emit: binned_fastas
    tuple val(meta), path("contig_bins.tsv"),  optional: true, emit: cluster_table
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    path "versions.yml"           , optional: true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args3 = task.ext.args3 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch contig_bins.tsv

    SemiBin2 \\
        $args \\
        single_easy_bin \\
        -i $fasta \\
        -b $bam \\
        -o $prefix \\
        -q \\
        -t $task.cpus \\
        $args3 >out 2>&1 || echo ""

    mv $prefix/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$( SemiBin2 --version )
    END_VERSIONS

"""
}
