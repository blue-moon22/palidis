/*
 * Convert FASTQ.GZ to FASTA, replace spaces in names with underscores, add f1/f2 and Seq number
 */
process convertToFasta {

	input:
    tuple val(sample_id), path(read1), path(read2)

	output:
	tuple val(sample_id), path("${sample_id}_1.fasta"), path("${sample_id}_2.fasta"), path("chunk*"), emit: fasta_ch
    tuple val(sample_id), path("${sample_id}_1.fasta"), path("${sample_id}_2.fasta"), emit: fasta_for_itr_ch

    script:
    fastq1 = "${sample_id}_1.fastq"
    fastq2 = "${sample_id}_2.fastq"

	"""
    for file in ${read1}
    do
        gunzip -c \$file >> ${fastq1}
    done

    for file in ${read2}
    do
        gunzip -c \$file >> ${fastq2}
    done

    convert_fastq_to_fasta.py -f ${fastq1} -r 1
    convert_fastq_to_fasta.py -f ${fastq2} -r 2
    cat ${sample_id}_*.fasta > ${sample_id}.fasta
    split -l 2000000 ${sample_id}.fasta chunk

    rm ${sample_id}.fasta
    rm ${fastq1} ${fastq2}
    """
}
