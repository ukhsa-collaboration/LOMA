#!/usr/bin/env python

# Calculates assembly statistics for each assembled genome in a directory.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
import csv
import sys


###### Functions

# Scaffold stats
def calculate_scaffold_stats(fasta_file):
	records = list(SeqIO.parse(fasta_file, "fasta"))
	count = len(records)
	total_asm = sum(len(record.seq) for record in records)
	lengths = [len(record.seq) for record in records]
	lengths.sort(reverse=True)
	gc_count = sum(record.seq.upper().count("G") + record.seq.upper().count("C") for record in records)
	all_count = sum(record.seq.upper().count("G") + record.seq.upper().count("C") + record.seq.upper().count("T") + record.seq.upper().count("A") for record in records)
	gc_cont = round((gc_count/all_count)*100, 2)
	s_n50 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.5), None)
	s_n90 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_asm * 0.9), None)
	gap_sum = sum(record.seq.count("N") for record in records)
	gap_perc = round((gap_sum/total_asm)*100, 2)
	return count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc


# Contig stats
def calculate_contig_stats(fasta_file, gap_size):
	records = [record for record in SeqIO.parse(fasta_file, "fasta")]
	new_records = []
	for record in records:
		seq = str(record.seq)
		contigs = seq.split("N" * int(gap_size)) # Split on gaps and create new records for split contigs
		for contig in contigs:
			if len(contig) > 0 and "N" not in contig:
				new_record = record[:]
				new_record.seq = Seq(contig)
				new_records.append(new_record)
	c_count = len(new_records)
	gap_count = len(new_records) - len(records)
	total_len = sum(len(record.seq) for record in new_records)
	lengths = [len(record.seq) for record in new_records]
	lengths.sort(reverse=True)
	c_n50 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_len * 0.5), None)
	c_n90 = next((length for length in lengths if sum(lengths[:lengths.index(length) + 1]) >= total_len * 0.9), None)
	return c_count, c_n50, gap_count, c_n90


# Write results to csv file
def write_stats_to_file(fasta_name, count, s_n50, c_count, c_n50, gap_count, gap_sum, total_asm, s_n90, c_n90, gc_cont, gap_perc, output, csvfile):
	sampleid = fasta_name.split(".")[0]
	binid = fasta_name.split(".")[1]
	writer = csv.writer(csvfile)
	writer.writerow([sampleid, binid, total_asm, count, s_n50, s_n90, c_count, c_n50, c_n90, gc_cont, gap_count, gap_sum, gap_perc])


# Parse arguments from the command line.
def parse_args():
	description = 'Calculates assembly statistics for each assembled genome in a directory. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--fasta_dir', help='Directory containing FASTA files (required)')
	parser.add_argument('--gap', help="Minimum gap length to be considered a scaffold (optional) [2]", default=2)
	parser.add_argument('--output', help="Output file prefix (optional) ['sample']", default="sample")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


###### Main
def main():
	args = parse_args()
	with open(args.output+'.assembly_stats.csv', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(["sampleid","bin", "assembly_length_bp", "scaffold_count", "scaffold_N50_bp", "scaffold_N90_bp","contig_count", "contig_N50_bp", "contig_N90_bp", "GC_perc", "gaps_count", "gaps_sum_bp", "gaps_perc"])

		# Find and process FASTA files
		for file_name in os.listdir(args.fasta_dir):
			if file_name.endswith('.fasta') or file_name.endswith('.fa'):
				file_path = os.path.join(args.fasta_dir, file_name)
				count, s_n50, gap_sum, total_asm, s_n90, gc_cont, gap_perc = calculate_scaffold_stats(file_path)
				c_count, c_n50, gap_count, c_n90 = calculate_contig_stats(file_path, args.gap)
				write_stats_to_file(file_name, count, s_n50, c_count, c_n50, gap_count, gap_sum, total_asm, s_n90, c_n90, gc_cont, gap_perc, args.output, csvfile)

if __name__ == "__main__":
	main()
