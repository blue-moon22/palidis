process palmem {

    input:
	tuple val(sample_id), path(fasta1), path(fasta2), path(db)

	output:
    tuple val(sample_id), path("${sample_id}_${db}_IR_1.fasta"), path("${sample_id}_${db}_paired_to_IR_2.fasta"), emit: ir_1_ch
    tuple val(sample_id), path("${sample_id}_${db}_paired_to_IR_1.fasta"), path("${sample_id}_${db}_IR_2.fasta"), emit: ir_2_ch
    tuple val(sample_id), path("${sample_id}_${db}_IR.tab"), emit: tab_ch

    script:
    min_itr_length = params.min_itr_length
    max_itr_length = params.max_itr_length
    kmer_length = params.kmer_length
	"""
	pal-mem -f1 ${fasta1} -f2 ${fasta2} -d ${db} -t ${task.cpus} -l ${min_itr_length} -m ${max_itr_length} -k ${kmer_length} -o ${sample_id}_${db}
	"""
}
