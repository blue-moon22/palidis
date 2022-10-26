#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Functions to get candidate ITRs and contigs with candidate insertion sequences.
"""

import argparse, sys
import re
from collections import defaultdict

def get_largest_index(fasta_file):
    """
    Function to get the largest sequence ID
    """

    with open(fasta_file, "r") as f:
        for line in f:
            if line[0] == '>':
                index = int(line.split("_")[0].replace(">Seq", ""))
    return(index)


def get_array_size(fasta_file1, fasta_file2):
    """
    Function to get the largest number of sequences.
    """

    size = max(get_largest_index(fasta_file1), get_largest_index(fasta_file2))
    return(size)


def set_positions(positions, left, right):
    """
    Function to return positions in ascending order
    """

    if positions:
        if left < positions[0]:
            positions[0] = left
        if right > positions[1]:
            positions[1] = right
    else:
        positions = [left, right]
    return positions


def get_ir_offset(tab_file, size):
    """
    Functiont to get the offset length of the inverted repeats (IRs) and assign
    to sequence ID
    """

    alloc1 = [0]*size
    alloc2 = [0]*size
    with open(tab_file, "r") as tab:
        for line in tab:
            seq1 = line.split('\t')[0]
            index = int(seq1.split("_")[0].replace("Seq", ""))
            left = int(seq1.split('_RCoord')[0].split('LCoord_')[1])
            right = int(seq1.split('_RCoord_')[1])
            if '_f1_' in seq1:
                alloc1[index-1] = set_positions(alloc1[index-1], left, right)
            elif '_f2_' in seq1:
                alloc2[index-1] = set_positions(alloc2[index-1], left, right)

            seq2 = line.split('\t')[1].replace('\n', '')
            left = int(seq2.split('_RCoord')[0].split('LCoord_')[1])
            right = int(seq2.split('_RCoord_')[1])
            index = int(seq2.split("_")[0].replace("Seq", ""))
            if '_f1_' in seq2:
                alloc1[index-1] = set_positions(alloc1[index-1], left, right)
            elif '_f2_' in seq2:
                alloc2[index-1] = set_positions(alloc2[index-1], left, right)
    return([alloc1, alloc2])


def write_contig_file(contigs, contigs_fasta_file, out_prefix):
    """
    Function to write contigs that contain candidate insertion sequences
    """

    # Write fasta file
    with open(out_prefix + '_contigs_with_candidate_itrs.fa', "w") as out_fasta:
        flag = 0
        with open(contigs_fasta_file, "r") as f:
            for line in f:
                if line[0] == ">":
                    contig = line.split("\n")[0].split(" ")[0].replace(">", "")
                    if contig in contigs:
                        out_fasta.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta.write(line)


def write_fasta_files(read_fasta1, read_fasta2, ir_arrays, position_arrays, out_prefix):
    """
    Function to write FASTA files that contain candidate ITRs.
    """

    # Get read names of IRs from first and second paired files
    # Write fasta file
    read_pairs_dict = {}
    with open(out_prefix + '_reads_with_candidate_itrs_1.fasta', "w") as out_fasta1:
        flag = 0
        with open(read_fasta1, "r") as f1:
            for line in f1:
                if line[0] == ">":
                    index = int(line.split("_")[0].replace(">Seq", ""))
                    if ir_arrays[0][index-1]:
                        if position_arrays[0][index-1]:
                            out_fasta1.write(line.replace('\n', '') + '_LCoord_' + str(position_arrays[0][index-1][0]) + '_RCoord_' + str(position_arrays[0][index-1][1]) + '\n')
                        else:
                            out_fasta1.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta1.write(line)

    with open(out_prefix + '_reads_with_candidate_itrs_2.fasta', "w") as out_fasta2:
        flag = 0
        with open(read_fasta2, "r") as f2:
            for line in f2:
                if line[0] == ">":
                    index = int(line.split("_")[0].replace(">Seq", ""))
                    if ir_arrays[1][index-1]:
                        if position_arrays[1][index-1]:
                            out_fasta2.write(line.replace('\n', '') + '_LCoord_' + str(position_arrays[1][index-1][0]) + '_RCoord_' + str(position_arrays[1][index-1][1]) + '\n')
                        else:
                            out_fasta2.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta2.write(line)


def split_match_flag(match_flag):
    """
    Function to find the exact positions of alignments from 'M' flags in sam
    files
    """
    split_matches = [0,0,0]
    cigar_dict = defaultdict(lambda: 0)
    if not any(cig in ['N', 'H', 'P', '=', 'X'] for cig in match_flag):
        codes = list(filter(None, re.split('[0-9]*', match_flag)))
        lengths = list(filter(None, re.split('M|S|I|D', match_flag)))
        count = 0
        for ind, c in enumerate(codes):
            if c == 'M':
                value = cigar_dict['M']
                cigar_dict['M'] = value + int(lengths[ind])
                count = 1
            elif c == 'S':
                if count:
                    split_matches[2] = int(lengths[ind])
                else:
                    split_matches[0] = int(lengths[ind])
                count = 1
            elif c == 'I':
                value = cigar_dict['M']
                cigar_dict['M'] = value - int(lengths[ind])
            elif c == 'D':
                value = cigar_dict['M']
                cigar_dict['M'] = value + int(lengths[ind])
        split_matches[1] = int(cigar_dict['M'])

    return split_matches


def process_sam_file(ir_positions_dict, sam_file, read_order, alloc_array):
    """
    Function to read information from the sam file and assign positions of IRs
    """

    with open(sam_file, "r") as sam:
        first_contig = ''

        for line in sam:
            current_contig = line.split("\t")[2]
            sam_flag = int(line.split('\t')[1])
            ir_paired_read = line.split("\t")[0]
            read_pair_id = ir_paired_read.split("_LCoord")[0][-3:]
            ir_position = 0
            index = int(ir_paired_read.split("_")[0].replace("Seq", ""))
            read_position = int(line.split("\t")[3])
            match_flag = line.split("\t")[5]
            split_matches = split_match_flag(match_flag)
            ir_offset = alloc_array[index-1]

            if read_pair_id == read_order[0]:
                if sam_flag in (99, 163, 97, 161, 65, 129, 67, 131):
                    start = read_position - split_matches[0] + ir_offset[0] - 1
                    ir_position = [(start, start + (ir_offset[1] - ir_offset[0])), True]
                elif sam_flag in (83, 147, 81, 145, 113, 177, 115, 179):
                    start = read_position - split_matches[0] + len(line.split("\t")[9]) - ir_offset[1]
                    ir_position = [(start, start + (ir_offset[1] - ir_offset[0])), True]

            # If it's the paired read e.g. Seq4_f2 of Seq_f1 and it's pair is unmapped (e.g. Seq_f1)
            elif read_pair_id == read_order[1]:
                if sam_flag in (73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137):
                    ir_position = [(read_position, read_position + len(line.split("\t")[9]) - 1), False] # Take the read position - not exact

            if ir_position:
                if current_contig in ir_positions_dict:
                    ir_positions_dict[current_contig][ir_paired_read] = ir_position
                else:
                    ir_positions_dict[current_contig] = {
                        ir_paired_read: ir_position
                    }

    return ir_positions_dict


def get_ir_positions(sam_files, alloc_arrays):
    """
    Function to get positions of IRs for both same files
    """

    # Get contigs and ITR cluster positions
    ir_positions_dict = {}
    ir_positions_dict = process_sam_file(ir_positions_dict, sam_files[0], ['_f1', '_f2'], alloc_arrays[0])
    ir_positions_dict = process_sam_file(ir_positions_dict, sam_files[1], ['_f2', '_f1'], alloc_arrays[1])

    return ir_positions_dict


def get_contigs_with_itrs(ir_positions_dict, position_arrays, size, MIN_IS_LEN, MAX_IS_LEN, out_prefix):
    """
    Function to write information on contigs that have IRs
    """

    alloc1 = [0]*size
    alloc2 = [0]*size
    contigs_w_ir = []
    sample_id = out_prefix.split('/')[len(out_prefix.split('/'))-1]
    with open(out_prefix + '_contigs_reads_ir_position_info.tab', "w") as out_tab:
        out_tab.write('sample_id\tcontig\tread\tIR\tstart_position\tend_position\n')

        for contig, read_info in ir_positions_dict.items():
            if len(read_info) > 1:
                reads = []
                start_positions = []
                end_positions = []
                ir = []
                for read, positions in read_info.items():
                    reads.append(read)
                    start_positions.append(positions[0][0])
                    end_positions.append(positions[0][1])
                    ir.append(positions[1])
                inds = [i for i in range(len(start_positions)) if any(abs(start_positions[i] - j) >= MIN_IS_LEN and abs(start_positions[i] - j) <= MAX_IS_LEN for j in start_positions)]
                for ind in inds:
                    index = int(reads[ind].split("_")[0].replace("Seq", ""))
                    if ir[ind]:
                        if '_f1' in reads[ind]:
                            alloc1[index-1] = 1
                            read_coord = '_LCoord_' + str(position_arrays[0][index-1][0]) + '_RCoord_' + str(position_arrays[0][index-1][1])
                        elif '_f2' in reads[ind]:
                            alloc2[index-1] = 1
                            read_coord = '_LCoord_' + str(position_arrays[1][index-1][0]) + '_RCoord_' + str(position_arrays[1][index-1][1])
                        out_tab.write(sample_id + '\t' + contig + '\t' + reads[ind] + read_coord + '\t' + str(ir[ind]) + '\t' + str(start_positions[ind]) + '\t' + str(end_positions[ind]) + '\n')
                    else:
                        read_coord = ''
                        if '_f1' in reads[ind]:
                            alloc2[index-1] = 1
                        elif '_f2' in reads[ind]:
                            alloc1[index-1] = 1

                    if contig not in contigs_w_ir:
                        contigs_w_ir.append(contig)

    return [alloc1, alloc2], contigs_w_ir


def get_arguments():

    parser = argparse.ArgumentParser(description='Get reads that contain candidate ITRs and contigs that associate with these.')
    parser.add_argument('--contig_fasta', '-c', dest='contig_fasta', required=True,
                        help='Input contig FASTA file.', type = str)
    parser.add_argument('--sam_file1', '-s1', dest='sam_file1', required=True,
                        help='Input sam file from mapping of reads to contigs.', type = str)
    parser.add_argument('--sam_file2', '-s2', dest='sam_file2', required=True,
                        help='Input sam file from mapping of reads to contigs.', type = str)
    parser.add_argument('--fasta1', '-f1', dest='fasta_file1', required=True,
                        help='Input first paired FASTA file.', type = str)
    parser.add_argument('--fasta2', '-f2', dest='fasta_file2', required=True,
                        help='Input second paired FASTA file.', type = str)
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input tab file of IR read pairs.', type = str)
    parser.add_argument('--min_is_len', '-min', dest='min_is_len', required=True,
                        help='Minimum length of insertion sequence', type = int, default = 500)
    parser.add_argument('--max_is_len', '-max', dest='max_is_len', required=True,
                        help='Maximum length of insertion sequence', type = int, default = 3000)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files.', type = str)
    return parser


def main(args):

    # Set min and max IS length
    MAX_IS_LEN = args.max_is_len
    MIN_IS_LEN = args.min_is_len

    # Allocate for fasta indexing
    size = get_array_size(args.fasta_file1, args.fasta_file2)

    print("Allocating empty arrays with IR size offset...")
    position_arrays = get_ir_offset(args.tab_file, size)

    # Get contigs, reads paired to reads with candidate ITRs and mapping positions
    print("Getting ITR positions...")
    contigs_read_and_positions = get_ir_positions([args.sam_file1, args.sam_file2], position_arrays)

    # Write tab file
    print("Writing tab file...")
    read_arrays, contigs_w_ir = get_contigs_with_itrs(contigs_read_and_positions, position_arrays, size, args.min_is_len, args.max_is_len, args.output_prefix)

    # Write contig FASTA file
    print("Writing contig FASTA with candidate ITRs...")
    write_contig_file(contigs_w_ir, args.contig_fasta, args.output_prefix)

    # Write FASTA files and get reads containing candidate ITRs
    print("Writing reads FASTAs with candidate ITRs...")
    write_fasta_files(args.fasta_file1, args.fasta_file2, read_arrays, position_arrays, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
