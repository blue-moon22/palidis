/*
 * Run Interproscan
 */
process runInterproscan {

    input:
    tuple val(sample_id), path(faa), path(db)

    output:
    tuple val(sample_id), path(faa), path("${faa}.tsv")

    script:
    """
    # Remove * from protein prediction
    sed -i 's/*//' ${faa}

    ./${db}/interproscan.sh -i ${faa} -f tsv -dp -cpu ${task.cpus}
    """
}
