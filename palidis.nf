/*
 * Nextflow pipeline for identifying insertion sequences from metagenomic data
 *
 * Author:
 * Victoria Carr vc11@sanger.ac.uk
 *
 */

nextflow.enable.dsl=2

// import modules
include { convertToFasta } from './modules/convertToFasta.nf'
include { filterContigs } from './modules/filterContigs.nf'
include { buildDB } from './modules/buildDB.nf'
include { palmem } from './modules/palmem.nf'
include { mapReads as mapReads1 } from './modules/mapreads.nf'
include { mapReads as mapReads2 } from './modules/mapreads.nf'
include { getCandidateITRs } from './modules/getCandidateITRs.nf'
include { clusterReads } from './modules/clusterReads.nf'
include { getITRs } from './modules/getITRs.nf'
include { runProdigal } from './modules/runProdigal.nf'
include { installInterproscan } from './modules/installInterproscan.nf'
include { runInterproscan } from './modules/runInterproscan.nf'
include { getISInfo } from './modules/getISInfo.nf'

workflow palidis {
    take:
    read_pair_ch
    contig_file_ch

    main:
    convertToFasta(read_pair_ch)

    palmem(convertToFasta.out)

    filterContigs(contig_file_ch)

    buildDB(filterContigs.out)

    /*
     * Map reads to contigs
     */
    Channel
    .of('1')
    .set { pair1_ch }

    palmem.out.ir_1_ch
    .join(buildDB.out.contig_db1_ch)
    .combine(pair1_ch)
    .set { contigs_reads1_ch }

    mapReads1(contigs_reads1_ch)

    Channel
    .of('2')
    .set { pair2_ch }

    palmem.out.ir_2_ch
    .join(buildDB.out.contig_db2_ch)
    .combine(pair2_ch)
    .set { contigs_reads2_ch }

    mapReads2(contigs_reads2_ch)

    /*
     * Get contigs and reads with candidate ITRs
     */
     contig_file_ch
     .join(mapReads1.out)
     .join(mapReads2.out)
     .join(convertToFasta.out)
     .join(palmem.out.tab_ch)
     .set { mapping_contigs_ch }

    getCandidateITRs(mapping_contigs_ch)

    clusterReads(getCandidateITRs.out.reads_itrs_ch)
    cluster_ch = clusterReads.out.cluster_ch

    cluster_ch
    .join(getCandidateITRs.out.tab_ch)
    .join(filterContigs.out)
    .set { into_get_itr_ch }

    getITRs(into_get_itr_ch)

    runProdigal(getITRs.out.is_fasta_for_prodigal_ch)
    installInterproscan()

    runProdigal.out
    .combine(installInterproscan.out)
    .set { proteins_ch }

    runInterproscan(proteins_ch)

    getITRs.out.is_tab_ch
    .join(runInterproscan.out)
    .join(getITRs.out.is_candidate_fasta_ch)
    .set { is_annot_ch }

    getISInfo(is_annot_ch)
    is_info_ch = getISInfo.out.txt
    is_fasta_ch = getISInfo.out.fasta

    emit:
    is_fasta_ch
    is_info_ch
}

workflow {
    // Define parameters
    batch_path = file("./${params.batch_name}")
    batch_path.mkdir()

    /*
     * Parameters
     */
    Channel
    .fromPath(params.manifest, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .map { row -> tuple(row.sample_id, file(row.read1), file(row.read2)) }
    .groupTuple()
    .set { read_pair_ch }

    Channel
    .fromPath(params.manifest, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .map { row -> tuple(row.sample_id, file(row.contigs_path)) }
    .set { contig_file_ch }

    palidis(read_pair_ch, contig_file_ch)

    // Publish IS fasta sequences
    palidis.out.is_fasta_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }

    // Publish annotations
    palidis.out.is_info_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }
}
