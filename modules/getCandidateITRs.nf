process getCandidateITRs {

    input:
    tuple val(sample_id), path(contig_file), path(sam_file1), path(sam_file2), path(fasta_ir1), path(fasta_ir2), path(tab_files)
    val c1
	val c2

	when:
	c1 == c2
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs_reads_ir_position_info.tab"), emit: tab_ch
    tuple val(sample_id), path("${sample_id}_reads_with_candidate_itrs.fasta"), emit: reads_itrs_ch

    script:
    min_is_len = params.min_is_len
    max_is_len = params.max_is_len

    """
    cat ${tab_files} > ${sample_id}.tab
    get_candidate_ITR_reads_and_IS_contigs.py \
        --contig_fasta ${contig_file} \
        --sam_file1 ${sam_file1} \
        --sam_file2 ${sam_file2} \
        --fasta1 ${fasta_ir1} \
        --fasta2 ${fasta_ir2} \
        --tab_file ${sample_id}.tab \
        --min_is_len ${min_is_len} \
        --max_is_len ${max_is_len} \
        --output_prefix ${sample_id}
    cat ${sample_id}_reads_with_candidate_itrs_1.fasta ${sample_id}_reads_with_candidate_itrs_2.fasta > ${sample_id}_reads_with_candidate_itrs.fasta
    rm ${sample_id}_reads_with_candidate_itrs_1.fasta ${sample_id}_reads_with_candidate_itrs_2.fasta
    """
}
