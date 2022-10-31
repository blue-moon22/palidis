#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Function to change headings in FASTA file
"""

import argparse, sys, os

def write_fasta(fasta_file, read):
    """
    This function writes a FASTA file from a FASTQ file and adds information to
    the headers
    """

    base_name = os.path.basename(fasta_file).split(".")[0]
    reversed = base_name[::-1]
    find = "_" + str(read)
    sample_name = reversed.replace(find[::-1], "", 1)[::-1]
    out_fasta_file = sample_name + "_" + str(read) + ".fasta"
    index = 1
    with open(out_fasta_file, "w") as out:
        with open(fasta_file, "r") as fasta:
            for line in fasta:
                if line[0] == '>':
                    new_id = ">Seq" + str(index) + "_f" + str(read)
                    out.write(new_id + "\n")
                    index += 1
                else:
                    out.write(line)


def get_arguments():
    parser = argparse.ArgumentParser(description='Convert FASTQ to FASTA with appropriate headers.')
    parser.add_argument('--fasta', '-f', dest='fasta', required=True,
                        help='Input FASTA file.', type = str)
    parser.add_argument('--read', '-r', dest='read', required=True,
                        help='Read file number (1 or 2).', type = str)
    return parser


def main(args):
    write_fasta(args.fasta, args.read)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
