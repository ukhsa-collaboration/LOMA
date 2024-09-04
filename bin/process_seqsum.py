#!/usr/bin/env python

# Fetches the sequencing summary file from Guppy output directories for parsing by pycoQC.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'

###### Imports

import os
import argparse


###### Functions

# Identify folders containing sequencing summary file and create a symlink
def create_sequencing_summary_symlink(input_dir, run_id):
	search_dir = os.path.join(input_dir, run_id)
	for root, _, files in os.walk(search_dir):
		for filename in files:
			if 'sequencing_summary' in filename:
				filepath = os.path.join(root, filename)
				try:
					os.symlink(filepath, os.path.basename(filepath))
				except OSError as e:
					print(f"Warning: Unable to create symlink for {filepath}: {e}")


# Parse arguments from the command line.
def parse_args():
	description = 'Concatenate identify sequencing summary files and create a symlink for processing by pycoQC. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("--input-dir", help="Input directory containing fastq files", required=True)
	parser.add_argument("--run_id", help="Run ID", default=None)
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


###### Main
def main():
	args = parse_args()
	create_sequencing_summary_symlink(args.input_dir, args.run_id)

if __name__ == "__main__":
	main()
