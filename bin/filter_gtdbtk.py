#!/usr/bin/env python

# Script to extract and relevant results from binned metagenome assembled genomes classified by GTDB-Tk (classify_wf). Designed to be run as part of a dedicated Nextflow pipeline.
# Relies on a custom GTDB definition table defined in the documentation

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports
import pandas as pd
import argparse
import os


###### Functions

def filter_gtdbtk(input_file, output_file, ANI_cutoff, align_fraction_cutoff, gtdb_definition_table, mode): # Process and subset GTDB-Tk results
	gtdbtk_results = pd.read_csv(input_file, sep='\t')
	gtdb_defs = pd.read_csv(gtdb_definition_table, sep='\t', dtype='unicode') # Import GTDB definition table

	gtdbtk_results = gtdbtk_results[(gtdbtk_results['closest_genome_ani']/100 >= ANI_cutoff) & (gtdbtk_results['closest_genome_af'] >= align_fraction_cutoff)] # Filter based on the specified cut-off values
	gtdbtk_results = gtdbtk_results[['user_genome','closest_genome_ani','closest_genome_af','classification']]
	merged_gtdbtk = gtdbtk_results.merge(gtdb_defs, left_on='classification', right_on='Original_ID').fillna('').astype(str)
	merged_gtdbtk = merged_gtdbtk[merged_gtdbtk["Target"] == 'Y'] # Subset only target species/assemblies
	merged_gtdbtk['fixed_id'] = merged_gtdbtk[['Clean_ID']].astype(str).replace(' ', '_',regex=True)

	if mode == "Nanopore":
		merged_gtdbtk['outstring'] = merged_gtdbtk[['user_genome', 'fixed_id','AMRFINDER','RESFINDER','MLST','KROCUS','gene_DB']].apply(lambda x: '.'.join(x), axis=1).replace(' ', '~', regex=True) + '.fasta' # Define output file name (to pass relevant metadata to channel)
	if mode == "Illumina":
		merged_gtdbtk['outstring'] = merged_gtdbtk[['user_genome', 'fixed_id','AMRFINDER','RESFINDER','MLST','SRST2','gene_DB']].apply(lambda x: '.'.join(x), axis=1).replace(' ', '~', regex=True) + '.fasta' # Define output file name (to pass relevant metadata to channel)
	print(merged_gtdbtk)
	return(merged_gtdbtk)

def extract_gtdbtk(merged_gtdbtk, output):
	for index, row in merged_gtdbtk.iterrows(): # Defined and create unique symlinks with renamed file names
		query_file = row['user_genome']
		input_name = f"{query_file}.fasta"
		symlink_name = row['outstring']
		wd = os.getcwd()
		input_path = os.path.join(wd,"input_bins", input_name)

		if not os.path.exists(symlink_name):
			if os.path.isfile(input_path):
				os.symlink(input_path,symlink_name)
			else:
				print("Input file missing, check fasta directory")

		else:
			print(f"Skipping symlink '{symlink_name}' as it already exists")

	merged_gtdbtk.to_csv(output+'.gtdb_filtered_report.txt', sep='\t', index=False) # Write to file


# Parse arguments from the command line.
def parse_args():
	description = 'Filter and extract data from a TSV file. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--input', required=True, help="Input TSV file")
	parser.add_argument('--output', required=True, help="Output TSV file")
	parser.add_argument('--ANI_cutoff', type=float, required=True, help="ANI cutoff value")
	parser.add_argument('--align_fraction_cutoff', type=float, required=True, help="ANI align fraction cutoff")
	parser.add_argument('--gtdb_definition_table', required=True, help="GTDB table")
	parser.add_argument('--mode', required=True, help="Sequencing type (Illumina or Nanopore) determines whether Krocus or SRST2 mode is included downstream")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()

###### Main
def main():
	args = parse_args()
	merged_gtdbtk = filter_gtdbtk(args.input, args.output, args.ANI_cutoff, args.align_fraction_cutoff, args.gtdb_definition_table, args.mode)
	extract_gtdbtk(merged_gtdbtk, args.output)

if __name__ == "__main__":
	main()
