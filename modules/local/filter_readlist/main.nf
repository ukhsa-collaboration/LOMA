process FILTER_READLIST {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fastq), path(minimap2_readlist, stageAs: "rl.1.txt"), path(kraken2_readlist, stageAs: "rl.2.txt")

    output:
    tuple val(meta), path(fastq), path("*.host_reads.txt"), optional: true, emit: host_readlist
    tuple val(meta), path(fastq), path("*.host_reads.original_IDs.txt"), optional: true, emit: host_readlist_originalIDs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat - $minimap2_readlist | cat - <(awk '\$1=="C"' $kraken2_readlist | cut -f2) | sort | uniq > ${prefix}.host_reads.txt
    cat ${prefix}.host_reads.txt | cut -f1 -d "_" | uniq > ${prefix}.host_reads.original_IDs.txt
    """
}
