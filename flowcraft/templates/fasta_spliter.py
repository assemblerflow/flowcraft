#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to split all fastas in a multifasta file into different
fasta files.

Code documentation
------------------

"""

import os
import sys


def main():

    cwd = os.getcwd()
    # a var to check if out_handle is started and if so it enables to control
    # how it should be closed
    out_handle = False
    # opens the input file of the process
    input_file = open(sys.argv[1])
    # a file with the list of all paths to fasta files that will be used by
    # fastANI
    list_files = open("files_fastani.txt", "w")
    # iterates by each entry in the fasta file
    for line in input_file:
        if line.startswith(">"):
            if out_handle:
                out_handle.close()
            # writes the output to fasta store folder inside cwd, respective
            # workdir
            path_to_file = os.path.join(cwd, "fasta_store",
                                        "_".join(line.split("_")[0:3])
                                        .replace(">", "") + ".fas")
            # writes to list of files
            list_files.write(path_to_file + "\n")
            out_handle = open(path_to_file, "w")
            out_handle.write(line)
        else:
            out_handle.write(line)

    out_handle.close()
    input_file.close()
    list_files.close()


if __name__ == "__main__":
    main()
