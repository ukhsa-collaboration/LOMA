#!/usr/bin/env python

# Fetches and concatenates basecalled FASTQ files from Guppy output directories.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'

###### Imports

import os
import argparse
import subprocess


###### Functions

# Identify relevant folders containing FASTQ files and concatenate them into a single *.fastq.gz output
def concatenate_fastq_files(input_dir, run_id, barcode, output_prefix, folder_name):
	output_file = output_prefix + ".concatenated.fastq.gz"
	print(output_file)
	search_dir = os.path.join(input_dir, run_id)
	keep_string = "fastq_pass"
	outputaa = set()

	with open(output_file, 'wb') as out_f:
		for root, dirs, files in os.walk(search_dir):
			if folder_name.lower() in root.lower():
				if keep_string.lower() in root.lower():
					if barcode is not None:
						barcode_string = "barcode" + barcode[-2:]
						if barcode_string.lower() not in root.lower():
							continue
						outputaa.add(root)
						if len(outputaa) > 1:
							raise ValueError("Multiple folders found, please ensure search strings and folder structures direct to a single *.fastq(.gz) directory")
						for filename in files:
							if filename.endswith('.fastq') or filename.endswith('.fastq.gz'):
								filepath = os.path.join(root, filename)
								if filename.endswith('.fastq.gz'):
									subprocess.run(['cat', filepath], stdout=out_f, check=True)

					elif barcode is None:
						outputaa.add(root)
						if len(outputaa) > 1:
							raise ValueError("Multiple folders found, please ensure search strings and folder structures direct to a single *.fastq(.gz) directory")
						for filename in files:
							if filename.endswith('.fastq') or filename.endswith('.fastq.gz'):
								filepath = os.path.join(root, filename)
								if 'barcode' in filepath:
									raise ValueError("Folder contains barcoded samples")
								if filename.endswith('.fastq.gz'):
									subprocess.run(['cat', filepath], stdout=out_f, check=True)


# Parse arguments from the command line.
def parse_args():
	description = 'Concatenate fastq files and create symlinks. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("--input-dir", help="Input directory containing fastq files", required=True)
	parser.add_argument("--run_id", help="Run ID", default=None)
	parser.add_argument("--barcode", help="Barcode string", default=None)
	parser.add_argument("--folder_name", help="basecaller (e.g. guppy, dorado etc)", default="guppy")
	parser.add_argument("--output_prefix", help="Output file prefix", required=True)
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


###### Main

def main():
	args = parse_args()
	concatenate_fastq_files(args.input_dir, args.run_id, args.barcode, args.output_prefix, args.folder_name)

if __name__ == "__main__":
		main()
