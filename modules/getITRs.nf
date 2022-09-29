/*
 * Get ITR information
 */
process getITRs {
    input:
    tuple val(sample_id), path(read_clstr_file), path(info_tab_file), path(assemblies_file)

    output:
    tuple val(sample_id), path("${sample_id}_insertion_sequence_annotations.tab"), emit: is_tab_ch
    tuple val(sample_id), path("${sample_id}_candidate_insertion_sequences.fasta"), emit: is_fasta_for_prodigal_ch
    tuple val(sample_id), path("${sample_id}_candidate_insertion_sequences.fasta"), emit: is_candidate_fasta_ch

    script:
    min_is_len = params.min_is_len
    max_is_len = params.max_is_len
    min_itr_len = params.min_itr_length
    max_itr_len = params.max_itr_length

    """
    assign_ITRs.py --cdhit_cluster_file ${read_clstr_file} \
        --info_tab_file ${info_tab_file} \
        --assemblies_fasta_file ${assemblies_file} \
        --min_is_len ${min_is_len} \
        --max_is_len ${max_is_len} \
        --min_itr_len ${min_itr_len} \
        --max_itr_len ${max_itr_len} \
        --cpus ${task.cpus} \
        --output_prefix ${sample_id}
    """
}
