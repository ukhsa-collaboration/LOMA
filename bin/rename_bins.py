#!/usr/bin/env python

# Rename FASTA files to a consistent format

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'

###### Imports

import os
import argparse

###### Functions

# Iterate through FASTA files, separate binned and unbinned MAGs, rename files (with bin ID or as 'unbinned') and contig number per-group
def rename_fasta_files_and_headers(input_dir, prefix):
	bin_count = 1 # Counter to ensure unique bin IDs

	for filename in os.listdir(input_dir):
		if filename.endswith(".fa") or filename.endswith(".fasta"): # Open FASTAs and create file name / header prefixes
			if filename.endswith('unbinned.fa'): # If unbinned
				new_filename = f"{prefix}.unbinned.fasta"
				header_prefix = f"{prefix}.unbinned.contig_"
			else: # If binned
				bx = ("%06d" % (bin_count,))
				new_filename = f"{prefix}.bin_{bx}.fasta"
				header_prefix = f"{prefix}.bin_{bx}.contig_"
				bin_count += 1
				contig_count = 1

			input_file = os.path.join(input_dir, filename)
			output_file = os.path.join(input_dir, new_filename)

			with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
				contig_count = 1 # Write contigs and headers to file incrementing contig number labels
				for line in infile:
					if line.startswith(">"):
						cx = ("%06d" % (contig_count,))
						outfile.write(f">{header_prefix}{cx}\n")
						contig_count += 1
					else:
						outfile.write(line)

			os.remove(input_file)

# Parse arguments from the command line.
def parse_args():
	description = 'Rename FASTA files to a consistent format. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("--input_dir",type=str, help="Input directory containing FASTA files")
	parser.add_argument("--prefix", type=str, help="Prefix for renaming files and headers")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


###### Main
def main():
	args = parse_args()
	rename_fasta_files_and_headers(args.input_dir, args.prefix)

if __name__ == "__main__":
	main()
