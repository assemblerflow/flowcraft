#!/usr/bin/env python3


import argparse
from itertools import groupby
import os


def replace_char(text):
    for ch in ['/', '`', '*', '{', '}', '[', ']', '(', ')', '#', '+', '-', '.', '!', '$', ':']:
        text = text.replace(ch, "_")
    return text

def getSequence(ref, fasta):

    entry = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

    for header in entry:
        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in entry.__next__())

        if ref == headerStr.replace('>',''):
            filename = os.path.join(os.getcwd(), ref.replace('/','_').split('|')[0])
            fasta_header = replace_char(headerStr)
            output_file = open(filename + '.fa', "w")
            output_file.write(">" + fasta_header + "\n" + seq.upper() + "\n")
            output_file.close()
            header_file = open("header.txt", "w")
            header_file.write(fasta_header)
            header_file.close()

def main():

    parser = argparse.ArgumentParser(prog='parse_fasta.py', description="Parse FASTA files for a specific header", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v0.1'))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-t', type=str, metavar='header of sequence to be retrieved',
                             help='Uncompressed fastq file containing mate 1 reads', required=True)
    parser_required.add_argument('-f', type=argparse.FileType('r'), metavar='/path/to/input/file.fasta',
                             help='Fasta with the sequences', required=True)

    args = parser.parse_args()

    getSequence(args.t, args.f)



if __name__ == "__main__":
    main()