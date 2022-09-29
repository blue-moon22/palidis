#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Functions to get the information for insertion sequences.
"""

import argparse, sys, os
import json
import re


def write_info(tab_file, prodigal_info, interpro_info, output_prefix):
    """
    Function to write the annotation information of the insertion sequences
    """

    is_name_dict = {}
    with open(f'{output_prefix}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\tinterpro_or_panther_accession\n")
        with open(tab_file, "r") as file:
            next(file)
            for line in file:
                is_name = line.split('\t')[0]

                # Interpro info
                if is_name in prodigal_info:
                    flag = 0
                    proteins = prodigal_info[is_name].keys()
                    for protein in proteins:
                        if protein in interpro_info:
                            accessions = []
                            for acc, items in interpro_info[protein].items():
                                length = is_name.split('_')[4]
                                transposase = '_'.join(re.split('-| ', items[0]))
                                transposase = ''.join(re.split('[^a-zA-Z0-9_]*', transposase))
                                protein_start = int(prodigal_info[is_name][protein][0])
                                start = str((protein_start - 1) + items[1][0])
                                end = str((protein_start - 1) + items[1][1])
                                if flag:
                                    new_is_name += f'_{transposase}_{start}-{end}'
                                else:
                                    new_is_name = f'IS_length_{length}_{transposase}_{start}-{end}'
                                    flag = 1
                                accessions.append(acc)

                    if flag:
                        out.write(f'{new_is_name}\t')
                        is_name_dict[is_name] = new_is_name
                        out.write('\t'.join(line.replace('\n', '').split('\t')[1:-1]) + '\t' + ';'.join(accessions) + '\n')

    return is_name_dict

def get_prodigal_info(aa_fasta):

    prodigal_dict = {}
    with open(aa_fasta) as faa:
        for line in faa:
            if line[0] == '>':
                header = line.split(' ')[0].replace('>', '')
                contig = header.rsplit('_', 1)[0]
                start = int(line.split(' ')[2])
                end = int(line.split(' ')[4])
                if contig in prodigal_dict:
                    prodigal_dict[contig][header] = (start, end)
                else:
                    prodigal_dict[contig] = {header: (start, end)}

    return prodigal_dict


def get_interproscan_info(interproscan_out):

    interpro_dict = {}
    with open(interproscan_out) as tsv:
        for line in tsv:
            if line.split('\t')[3] == 'PANTHER':
                annotation = line.split('\t')[5]
            else:
                annotation = line.split('\t')[12].replace('\n', '')

            if any(x in ''.join(re.split('[^a-zA-Z0-9]*', annotation)).lower() for x in ['transposase', 'integraselike', 'ribonucleaseh']):
                protein = line.split('\t')[0]
                start = (int(line.split('\t')[6])-1)*3
                end = (int(line.split('\t')[7])-1)*3

                if line.split('\t')[3] == 'PANTHER':
                    accession = line.split('\t')[4]
                else:
                    accession = line.split('\t')[11]

                if protein in interpro_dict:
                    interpro_dict[protein][accession] = [annotation, (start, end)]
                else:
                    interpro_dict[protein] = {accession: [annotation, (start, end)]}

    return interpro_dict


def write_fasta_file(fasta_file, is_name_dict, output_prefix):

    flag = 0
    with open(f'{output_prefix}_insertion_sequences.fasta', "w") as out:
        with open(fasta_file, 'r') as file:
            for line in file:
                if line[0] == '>':
                    is_name = line.replace('>', '').replace('\n', '')
                    if is_name in is_name_dict:
                        new_name = is_name_dict[is_name]
                        out.write(f'>{new_name}\n')
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out.write(line)


def get_arguments():

    parser = argparse.ArgumentParser(description='Get tab file of annotations.')
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input "insertion_sequence_annotations.tab" file.', type = str)
    parser.add_argument('--aa_fasta', '-p', dest='aa_fasta', required=True,
                        help='Input "_insertion_sequences.faa file."', type = str)
    parser.add_argument('--interproscan_out', '-i', dest='interproscan_out', required=True,
                        help='Input "_insertion_sequences.faa.tsv" file."', type = str)
    parser.add_argument('--fasta_file', '-f', dest='fasta_file', required=True,
                        help='FASTA file of candidate insertion sequences.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    prodigal_info_dict = get_prodigal_info(args.aa_fasta)
    interpro_info_dict = get_interproscan_info(args.interproscan_out)

    is_name_dict = write_info(args.tab_file, prodigal_info_dict, interpro_info_dict, args.output_prefix)
    write_fasta_file(args.fasta_file, is_name_dict, args.output_prefix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
