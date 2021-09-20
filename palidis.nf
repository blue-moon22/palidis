/*
 * Nextflow pipeline for identifying insertion sequences from metagenomic data
 *
 * Author:
 * Victoria Carr <victoriacarr018@gmail.com
 *
 */

nextflow.enable.dsl=2

/*
 * Convert FASTQ.GZ to FASTA, replace spaces in names with underscores, add f1/f2 and Seq number
 */
process convertToFasta {

	input:
    tuple val(sample_id), path(read1), path(read2)

	output:
	tuple val(sample_id), path("${sample_id}_1.fasta"), path("${sample_id}_2.fasta")

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

    rm ${fastq1} ${fastq2}
    """
}

/*
 * Run pal-MEM
 */
process palmem {

    input:
	tuple val(sample_id), path(fasta1), path(fasta2)

	output:
    tuple val(sample_id), path("${sample_id}_IR_1.fasta"), path("${sample_id}_paired_to_IR_2.fasta"), emit: ir_1_ch
    tuple val(sample_id), path("${sample_id}_paired_to_IR_1.fasta"), path("${sample_id}_IR_2.fasta"), emit: ir_2_ch

    script:
    min_itr_length = params.min_itr_length
    kmer_length = params.kmer_length
    split = params.split
	"""
	pal-mem -f1 ${fasta1} -f2 ${fasta2} -t ${task.cpus} -l ${min_itr_length} -k ${kmer_length} -o ${sample_id} -d ${split}
	"""
}

/*
 * Build bowtie2 database
 */
process buildDB {

    input:
    tuple val(sample_id), path(contigs_path)

    output:
    tuple val(sample_id), path(contigs_path), path("contigs_db.tar"), emit: contig_db1_ch
    tuple val(sample_id), path("contigs_db.tar"), emit: contig_db2_ch

    """
    bowtie2-build ${contigs_path} ${sample_id}_contigs
    tar -cf contigs_db.tar ${sample_id}_contigs*
    """
}

process mapReads1 {

	input:
	tuple val(sample_id), path(fasta1), path(fasta2), path(contigs_path), path(db_path)

	output:
	tuple val(sample_id), path(contigs_path), path(fasta1), path("${prefix}.sam.mapped.sorted")

    script:
    prefix = "${sample_id}_itr_in_1"
	"""
    tar -xf ${db_path}
	bowtie2 --very-sensitive-local -x ${sample_id}_contigs -1 ${fasta1} -2 ${fasta2} -S ${prefix}.sam -p ${task.cpus} -f
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

process mapReads2 {

	input:
	tuple val(sample_id), path(fasta1), path(fasta2), path(db_path)

	output:
	tuple val(sample_id), path(fasta2), path("${prefix}.sam.mapped.sorted")

    script:
    prefix = "${sample_id}_itr_in_2"
	"""
    tar -xf ${db_path}
	bowtie2 --very-sensitive-local -x ${sample_id}_contigs -1 ${fasta1} -2 ${fasta2} -S ${prefix}.sam -p ${task.cpus} -f
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

process getCandidateITRs {

    input:
    tuple val(sample_id), path(contig_file), path(fasta_ir1), path(sam_file1), path(fasta_ir2), path(sam_file2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs_reads_ir_position_info.tab"), emit: tab_ch
    tuple val(sample_id), path("${sample_id}_reads_with_candidate_itrs.fasta"), emit: reads_itrs_ch

    """
    cat ${sam_file1} ${sam_file2} > combined.sam
    get_insert_size.py combined.sam > sam.stats
    insert_size=\$(grep 'Read span' sam.stats | cut -f4 -d' ' | sed 's/,//')
    read_length=\$(grep 'Read length' sam.stats | cut -f4 -d' ' | sed 's/,//')

    get_candidate_ITR_reads_and_IS_contigs.py \
        --contig_fasta ${contig_file} \
        --sam_file1 ${sam_file1} \
        --sam_file2 ${sam_file2} \
        --fasta1 ${fasta_ir1} \
        --fasta2 ${fasta_ir2} \
        --insert_size \${insert_size} \
        --read_length \${read_length} \
        --output_prefix ${sample_id}
    cat ${sample_id}_reads_with_candidate_itrs_1.fasta ${sample_id}_reads_with_candidate_itrs_2.fasta > ${sample_id}_reads_with_candidate_itrs.fasta
    rm ${sample_id}_reads_with_candidate_itrs_1.fasta ${sample_id}_reads_with_candidate_itrs_2.fasta
    """
}

/*
 * Cluster reads
 */
process clusterReads {

    input:
    tuple val(sample_id), file(read_file)

    output:
    tuple val(sample_id), path("${output_prefix}.fasta.clstr"), path("${sample_id}_irs.fasta"), emit: nonred_read_ch

    script:
    G=params.cd_hit_G
    aL=params.cd_hit_aL
    aS=params.cd_hit_aS
    A=params.min_itr_length
    output_prefix="${sample_id}_nonred_G${G}_aL${aL}_aS${aS}_A${A}"

    """
    clip_reads.py --read_fasta ${read_file} --output_prefix ${sample_id}
    cd-hit-est -i ${sample_id}_irs.fasta -o ${output_prefix}.fasta -G ${G} -aL ${aL} -aS ${aS} -A ${A} -M 64000 -T ${task.cpus} -d 0
    """
}

/*
 * Get ITR information
 */
process getITRs {
    input:
    tuple val(sample_id), path(read_clstr_file), path(clipped_fasta), path(info_tab_file)

    output:
    tuple val(sample_id), path("${sample_id}_contigs_reads_itr_position_info.tab"), emit: itr_tab_ch
    path("${sample_id}_ITRs.fasta"), emit: itr_fasta_ch

    """
    assign_ITRs.py --clipped_reads ${clipped_fasta} --cdhit_cluster_file ${read_clstr_file} --info_tab_file ${info_tab_file} --output_prefix ${sample_id}
    """
}

/*
 * Collect IS annotations
 */
process collectAnnotations {
    input:
    tuple val(sample_id), path(itr_pos_info)

    output:
    path("${sample_id}_insertion_sequence_annotations.tab"), emit: is_tab_ch
    path("${sample_id}_reads_itr_clusters.txt"), emit: itr_clusters_ch

    script:
    collect_annotations_script = file(params.collect_annotations_script)

    """
    Rscript --vanilla ${collect_annotations_script} --input ${itr_pos_info} --output ${sample_id}
    """
}

/*
 * Combine ITR sequences
 */
process createITRCatalog {
    input:
    path itr_fastas

    output:
    path "all_ITRs.fasta", emit: all_itrs_ch

    script:
    """
    cat ${itr_fastas} > all_ITRs.fasta
    seqtk seq -r all_ITRs.fasta > all_ITRs_rev.fasta
    cat all_ITRs_rev.fasta >> all_ITRs.fasta
    """
}

/*
 * Run pal-MEM on ITRs
 */
process palmemOnITRs {

    input:
	path all_itrs

	output:
    path "all_IR.tab", emit: all_itrs_tab

    script:
    min_itr_length = params.min_itr_length
    kmer_length = params.kmer_length
    split = params.split
	"""
	pal-mem -fu ${all_itrs} -t ${task.cpus} -l ${min_itr_length} -k ${kmer_length} -o all -d ${split}
	"""
}

workflow get_IS_annotations {
    take:
    read_pair_ch
    contig_file_ch

    main:
    convertToFasta(read_pair_ch)
    fasta_ch = convertToFasta.out

    palmem(fasta_ch)

    buildDB(contig_file_ch)

    /*
     * Map reads to contigs
     */
    palmem.out.ir_1_ch
    .join(buildDB.out.contig_db1_ch)
    .set { contigs_reads1_ch }

    mapReads1(contigs_reads1_ch)
    mapping1_ch = mapReads1.out

    palmem.out.ir_2_ch
    .join(buildDB.out.contig_db2_ch)
    .set { contigs_reads2_ch }

    mapReads2(contigs_reads2_ch)
    mapping2_ch = mapReads2.out

    /*
     * Get contigs and reads with candidate ITRs
     */
     mapping1_ch
     .join(mapping2_ch)
     .set { mapping_contigs_ch }

    getCandidateITRs(mapping_contigs_ch)

    reads_itrs_ch = getCandidateITRs.out.reads_itrs_ch
    tab_ch = getCandidateITRs.out.tab_ch

    clusterReads(reads_itrs_ch)
    nonred_read_ch = clusterReads.out.nonred_read_ch

    nonred_read_ch
    .join(tab_ch)
    .set { into_get_itr_ch }

    getITRs(into_get_itr_ch)
    itr_tab_ch = getITRs.out.itr_tab_ch
    itr_fasta_ch = getITRs.out.itr_fasta_ch

    collectAnnotations(itr_tab_ch)
    itr_clusters_ch = collectAnnotations.out.itr_clusters_ch
    is_tab_ch = collectAnnotations.out.is_tab_ch

    emit:
    itr_fasta_ch
    itr_clusters_ch
    is_tab_ch
}

workflow cluster_ITRs {
    take:
    itrs_fasta_ch
    is_annot_ch

    main:
    createITRCatalog(itrs_fasta_ch)
    all_itrs_ch = createITRCatalog.out

    palmemOnITRs(all_itrs_ch)
    all_itrs_tab_ch = palmemOnITRs.out

    // Create cluster membership in R using https://stackoverflow.com/questions/45736832/simple-network-cluster-membership-from-two-column-data-frame and join to each dataframe
    emit:
    all_itrs_tab_ch
}

workflow {
    // Define parameters
    batch_path = file("./${params.batch_name}")
    batch_path.mkdir()

    if (params.get_IS_annotations) {
        /*
         * Parameters
         */
        Channel
        .fromPath(params.manifest)
        .splitCsv(header:true, sep:"\t")
        .map { row -> tuple(row.sample_id, file(row.read1), file(row.read2)) }
        .groupTuple()
        .set { read_pair_ch }

        Channel
        .fromPath(params.manifest)
        .splitCsv(header:true, sep:"\t")
        .map { row -> tuple(row.sample_id, file(row.contigs_path)) }
        .set { contig_file_ch }

        get_IS_annotations(read_pair_ch, contig_file_ch)

        // Publish batch of candidate ITRs
        get_IS_annotations.out.itr_fasta_ch
        .flatten()
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }

        // Publish itr clusters file for batch
        get_IS_annotations.out.itr_clusters_ch
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }

        // Publish tab file for batch
        get_IS_annotations.out.is_tab_ch
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }
    }

    if (params.cluster_ITRs) {

        Channel
        .fromPath("${batch_path}/*_ITRs.fasta", checkIfExists:true)
        .collect()
        .set { itrs_fasta_ch }

        Channel
        .fromPath("${batch_path}/*_insertion_sequence_annotations.tab", checkIfExists:true)
        .collect()
        .set { is_annot_ch }

        cluster_ITRs(itrs_fasta_ch, is_annot_ch)
    }
}
