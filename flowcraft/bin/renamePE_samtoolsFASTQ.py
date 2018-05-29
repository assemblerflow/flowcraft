#!/usr/bin/env python2

#TODO - change to py3
# -*- coding: utf-8 -*-

"""
renamePE_samtoolsFASTQ.py - Rename the fastq headers with PE terminations
that were not include in samtools fastq command
<https://github.com/miguelpmachado/pythonScripts>
Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>
Last modified: January 10, 2017
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import time
import argparse
import itertools


version = '0.1'


def formartFastqHeaders(in_fastq_1, in_fastq_2, outdir):
	out_fastq_1 = os.path.join(outdir, os.path.splitext(os.path.basename(in_fastq_1))[0] + '.headersRenamed_1.fq')
	out_fastq_2 = os.path.join(outdir, os.path.splitext(os.path.basename(in_fastq_2))[0] + '.headersRenamed_2.fq')
	writer_in_fastq_1 = open(out_fastq_1, 'wt')
	writer_in_fastq_2 = open(out_fastq_2, 'wt')
	outfiles = [out_fastq_1, out_fastq_2]
	with open(in_fastq_1, 'rtU') as reader_in_fastq_1, open(in_fastq_2, 'rtU') as reader_in_fastq_2:
		plus_line = True
		quality_line = True
		number_reads = 0
		for in_1, in_2 in itertools.izip(reader_in_fastq_1, reader_in_fastq_2):
			if len(in_1) > 0:
				in_1 = in_1.splitlines()[0]
				in_2 = in_2.splitlines()[0]
				if in_1.startswith('@') and plus_line and quality_line:
					if in_1 != in_2:
						sys.exit('The PE fastq files are not aligned properly!')
					in_1 += '/1' + '\n'
					in_2 += '/2' + '\n'
					writer_in_fastq_1.write(in_1)
					writer_in_fastq_2.write(in_2)
					plus_line = False
					quality_line = False
				elif in_1.startswith('+') and not plus_line:
					in_1 += '\n'
					writer_in_fastq_1.write(in_1)
					writer_in_fastq_2.write(in_1)
					plus_line = True
				elif plus_line and not quality_line:
					in_1 += '\n'
					in_2 += '\n'
					writer_in_fastq_1.write(in_1)
					writer_in_fastq_2.write(in_2)
					writer_in_fastq_1.flush()
					writer_in_fastq_2.flush()
					number_reads += 1
					quality_line = True
				else:
					in_1 += '\n'
					in_2 += '\n'
					writer_in_fastq_1.write(in_1)
					writer_in_fastq_2.write(in_2)
	return number_reads, outfiles


def compressionType(file_to_test):
	magic_dict = {'\x1f\x8b\x08': ['gzip', 'gunzip'], '\x42\x5a\x68': ['bzip2', 'bunzip2']}

	max_len = max(len(x) for x in magic_dict)

	with open(file_to_test, 'r') as reader:
		file_start = reader.read(max_len)

	for magic, filetype in magic_dict.items():
		if file_start.startswith(magic):
			return filetype
	return None


def runTime(start_time):
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print 'Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's'
	return time_taken


def main():
	parser = argparse.ArgumentParser(prog='renamePE_samtoolsFASTQ.py', description='Rename the fastq headers with PE terminations that were not include in samtools fastq command', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-1', '--fastq_1', type=argparse.FileType('r'), metavar='/path/to/input/file_1.fq', help='Uncompressed fastq file containing mate 1 reads', required=True)
	parser_required.add_argument('-2', '--fastq_2', type=argparse.FileType('r'), metavar='/path/to/input/file_2.fq', help='Uncompressed fastq file containing mate 2 reads', required=True)

	parser_optional_general = parser.add_argument_group('General facultative options')
	parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')

	args = parser.parse_args()

	print '\n' + 'STARTING renamePE_samtoolsFASTQ.py' + '\n'
	start_time = time.time()

	fastq_files = [os.path.abspath(args.fastq_1.name), os.path.abspath(args.fastq_2.name)]

	print 'Check if files are compressed' + '\n'
	for fastq in fastq_files:
		if compressionType(fastq) is not None:
			sys.exit('Compressed fastq files found')

	outdir = os.path.abspath(args.outdir)
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	print 'Renaming fastq headers' + '\n'
	number_reads, outfiles = formartFastqHeaders(fastq_files[0], fastq_files[1], outdir)

	print 'It was written ' + str(number_reads) + ' read pairs in ' + str(outfiles) + ' files' + '\n'

	print '\n' + 'END renamePE_samtoolsFASTQ.py'
	time_taken = runTime(start_time)
	del time_taken


if __name__ == "__main__":
	main()