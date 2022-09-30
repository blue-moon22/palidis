/*
 * Get IS Information
 */
process getISInfo {

    input:
    tuple val(sample_id), path(is_tab_file), path(faa_file), path(tsv_file), path(fasta_file)

    output:
    path("${output_txt}"), emit: txt
    path("${output_fasta}"), emit: fasta

    script:
    output_txt="${sample_id}_insertion_sequences_info.txt"
    output_fasta="${sample_id}_insertion_sequences.fasta"
    """
    get_IS_info.py \
        --tab_file ${is_tab_file} \
        --aa_fasta ${faa_file} \
        --interproscan_out ${tsv_file} \
        --fasta_file ${fasta_file} \
        --output_prefix ${sample_id}
    """
}
