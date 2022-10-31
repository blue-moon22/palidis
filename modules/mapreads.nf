process mapReads {

	input:
	tuple val(sample_id), path(fastas1), path(fastas2), path(db_path), val(pair)
	val c1
	val c2

	when:
	c1 == c2

	output:
	tuple val(sample_id), path("${prefix}.sam.mapped.sorted")

    script:
    prefix = "${sample_id}_itr_in_${pair}"
	"""
	cat ${fastas1} > tmp_1.fasta
	sed -e '/^>/s/\$/@/' -e 's/^>/#/' tmp_1.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort | uniq | sed -e 's/^/>/' -e 's/\\t/\\n/' | sed '1d' > ${sample_id}_1.fasta
	cat ${fastas2} > tmp_2.fasta
	sed -e '/^>/s/\$/@/' -e 's/^>/#/' tmp_2.fasta | tr -d '\n' | tr "#" "\n" | tr "@" "\t" | sort | uniq | sed -e 's/^/>/' -e 's/\\t/\\n/' | sed '1d' > ${sample_id}_2.fasta

    tar -xf ${db_path}
	bowtie2 --very-sensitive-local -x ${sample_id}_contigs -1 ${sample_id}_1.fasta -2 ${sample_id}_2.fasta -S ${prefix}.sam -p ${task.cpus} -f
	samtools view -S -b ${prefix}.sam -@ ${task.cpus} > ${prefix}.bam
	rm ${prefix}.sam
	samtools view -b -F 4 ${prefix}.bam -@ ${task.cpus} > ${prefix}.bam.mapped
	rm ${prefix}.bam
	samtools sort ${prefix}.bam.mapped -o ${prefix}.bam.mapped.sorted -@ ${task.cpus}
	rm ${prefix}.bam.mapped
	samtools index ${prefix}.bam.mapped.sorted
	samtools view ${prefix}.bam.mapped.sorted > ${prefix}.sam.mapped.sorted
	rm ${prefix}.bam.mapped.sorted*
    rm ${sample_id}_contigs*
	"""
}
