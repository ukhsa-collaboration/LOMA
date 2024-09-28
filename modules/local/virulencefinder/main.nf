process VIRULENCEFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/virulencefinder:2.0.4--hdfd78af_0' :
        'https://depot.galaxyproject.org/singularity/virulencefinder%3A2.0.4--hdfd78af_0'  }"

    input:
    tuple val(meta), path(assembly)
    path(vfdb)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    virulencefinder.py -i $assembly -x -o ./ -p $vfdb
    mv data.json ${prefix}.virulencefinder.json
    mv results_tab.tsv ${prefix}.virulencefinder.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rename_bins.py: \$(rename_bins.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS

    """
}
