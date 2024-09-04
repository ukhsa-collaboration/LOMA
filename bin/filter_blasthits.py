#!/usr/bin/env python

# Script to summarize BLAST results

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports

import csv
import argparse
from collections import defaultdict


###### Functions

# Process the workflow files into a format that's easy to match with the blastn results
def process_workflow(workflow):
	workflow_data = {}
	with open(workflow, 'r', newline='') as workflow_input:
		reader = csv.reader(workflow_input, delimiter='\t')
		for row in reader:
			gene = row[0]
			ident = float(row[1])
			description = row[4]
			workflow_data[gene] = {'ident': ident, 'description': description}

	return(workflow_data)

# Filter the blastn results using the processed workflow data
def filter_results(blast_results, workflow_data, overlap_len):
	filtered_rows = []
	seen_hits = set()
	gene_count = defaultdict(int)
	with open(blast_results, 'r', newline='') as blast_input:
		reader = csv.reader(blast_input, delimiter='\t')
		for row in reader:
			sseqid = row[1]
			pident = float(row[2])
			qseqid = row[0]
			qstart = int(row[6])
			qend = int(row[7])
			length = int(row[3])
			qcovs = int(row[12])

			# Identify full length blastn hits above the specified pident
			if sseqid in workflow_data and pident >= workflow_data[sseqid]['ident'] and (length/qcovs)*100 == 100:
				# Check for overlapping hits and keep the first one
				hit_key = (qseqid, sseqid, qstart, qend)
				if not any(existing_hit for existing_hit in seen_hits if overlap(existing_hit, hit_key, overlap_len)):
					seen_hits.add(hit_key)
					row.append(workflow_data[sseqid]['description'])
					row.insert(14, "Good")
					filtered_rows.append(row)
					gene_count[sseqid] += 1

			# Identify partial/longer length blastn hits above the specified pident
			if sseqid in workflow_data and pident >= workflow_data[sseqid]['ident'] and (length/qcovs)*100 >= pident and (length/qcovs)*100 < 100:
				# Check for overlapping hits and keep the first one
				hit_key = (qseqid, sseqid, qstart, qend)
				if not any(existing_hit for existing_hit in seen_hits if overlap(existing_hit, hit_key, overlap_len)):
					seen_hits.add(hit_key)
					row.append(workflow_data[sseqid]['description'])
					row.insert(14, "Uncertain")
					filtered_rows.append(row)
					gene_count[sseqid] += 1

		return(filtered_rows, gene_count)


# Write filtered hits to files
def write_outputs(filtered_rows, gene_count, workflow_data, output):
	# Write filtered rows to output file
	output_csv_path = output + '.blastn_filtered_hits.tsv'
	with open(output_csv_path, 'w', newline='') as output_csv:
		writer = csv.writer(output_csv, delimiter='\t')
		writer.writerow(['qseqid','sseqid','pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore','slen','Description', 'Quality'])
		for row in filtered_rows:
			writer.writerow(row)

	# Write gene counts to output file
	output_count_csv_path = output + '.blastn_summary.tsv'
	with open(output_count_csv_path, 'w', newline='') as count_csv:
		count_writer = csv.writer(count_csv, delimiter='\t')
		count_writer.writerow(['gene', 'description','count'])
		for gene, count in gene_count.items():
			count_writer.writerow([gene, workflow_data[gene]['description'],count])

# Define how overlaps are evaluated
def overlap(hit1, hit2, overlap_len):
	return( hit1[0] == hit2[0] and hit1[1] == hit2[1] and max(0, min(hit1[3], hit2[3]) - max(hit1[2], hit2[2])) > overlap_len)


# Parse arguments from the command line.
def parse_args():
	description = 'Filter and summarize blastn results. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--blast_results', required=True, help="Blast results")
	parser.add_argument('--workflow', required=True, help="Workflow.txt file")
	parser.add_argument('--overlap', required=False, default=3, type=float, help="Minimum overlap (in bp) required to count as a single hit [3]")
	parser.add_argument('--output', required=True, help="Output TSV prefix")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))
	
	return parser.parse_args()


###### Main
def main():
	args = parse_args()

	workflow_data = process_workflow(args.workflow)
	filtered_rows, gene_count = filter_results(args.blast_results, workflow_data, args.overlap)
	write_outputs(filtered_rows, gene_count, workflow_data, args.output)


if __name__ == "__main__":
	main()
