process RESFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://genomicepidemiology/resfinder:4.5.0' :
        'biocontainers/resfinder:4.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.txt"), emit: report
    tuple val(meta), path("*.ResFinder_results.txt"), emit: results
    tuple val(meta), path("*.ResFinder_results_tab.txt"), emit: results_resfinder
    tuple val(meta), path("*.PointFinder_results.txt"), optional: true, emit: results_pointfinder

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    organism_param_1 = meta.containsKey("resfinder") ? "-s" : "-s"
    organism_param_2 = !(meta.resfinder.matches("")) ? "${meta.resfinder}" : "Other"
    if(params.RESFINDER.db) { resf_db = "export CGE_RESFINDER_RESGENE_DB=${params.RESFINDER.db}"} else {resf_db = ""}  
    if(params.RESFINDER.db) { pointf_db = "export CGE_RESFINDER_RESPOINT_DB=${params.POINTFINDER.db}"} else {pointf_db = ""}
    """
    $resf_db
    $pointf_db

    python \\
      -m resfinder \\
      -ifa $assembly \\
      -o ${prefix} \\
      --acquired \\
      --point ${organism_param_1} "${organism_param_2}" $args 

    mv */*.txt ./
    mv ResFinder_results_tab.txt ${prefix}.ResFinder_results_tab.txt
    mv ResFinder_results.txt ${prefix}.ResFinder_results.txt
    mv PointFinder_results.txt ${prefix}.PointFinder_results.txt 2>/dev/null || echo ""

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$(echo \$(python -m resfinder --version 2>&1))
    END_VERSIONS
    """
}
