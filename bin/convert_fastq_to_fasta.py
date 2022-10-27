#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Function to convert a FASTQ file into a FASTA file with additional information.
"""

import argparse, sys, os

def write_fasta(fastq_file, read):
    """
    This function writes a FASTA file from a FASTQ file and adds information to
    the headers
    """

    base_name = os.path.basename(fastq_file).split(".")[0]
    reversed = base_name[::-1]
    find = "_" + str(read)
    sample_name = reversed.replace(find[::-1], "", 1)[::-1]
    fasta_file = sample_name + "_" + str(read) + ".fasta"
    index = 1
    with open(fasta_file, "w") as out:
        with open(fastq_file, "r") as fastq:
            for line in fastq:
                if (index - 1) % 4 == 0:
                    new_id = ">Seq" + str(index) + "_nstart_" + sample_name + "_nend_" + line.replace('\n', "").replace(" ", "_") + "_f" + str(read)
                    out.write(new_id + "\n")
                elif (index - 1) % 4 == 1:
                    out.write(line)
                index += 1


def get_arguments():
    parser = argparse.ArgumentParser(description='Convert FASTQ to FASTA with appropriate headers.')
    parser.add_argument('--fastq', '-f', dest='fastq', required=True,
                        help='Input FASTQ file.', type = str)
    parser.add_argument('--read', '-r', dest='read', required=True,
                        help='Read file number (1 or 2).', type = str)
    return parser


def main(args):
    write_fasta(args.fastq, args.read)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
