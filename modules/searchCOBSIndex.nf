process searchCOBSIndex {

    container 'leandroishilima/661k_query_indexes:0.0.1'

    input:
    tuple val(sample_id), file(query), file(cobs_index)

    output:
    tuple val(sample_id), path("${output}"), optional: true

    script:
    cobs_threshold=params.cobs_threshold
    output="${query}_${cobs_threshold}_results_table.txt"

    """
    if [ -s ${query} ]
    then
        # query COBS index
        cobs query -i ${cobs_index} -f ${query} -t ${cobs_threshold} > ${query}_${cobs_threshold}_results.txt

        # run samtools faidx to get length of each query sequence
        samtools faidx ${query}

        # calculate percentage of kmers present rather than number of kmers present
        cobs_search_result_to_table.py --cobs_outfile ${query}_${cobs_threshold}_results.txt --fai_file ${query}.fai --outname ${output}
    fi
    """
}
