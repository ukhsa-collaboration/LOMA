#!/usr/bin/env python

# Merges and summarizes targeted gene finding results, BLASTN, GeneFinder and VirulenceFinder (optional)

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports

import pandas as pd
import xml.etree.ElementTree as ET
import argparse


###### Functions

# Process GeneFinder hits
def extract_genefinder_hits(genefinder_results):

	tree = ET.parse(genefinder_results)
	root = tree.getroot()

	gene_names = set()

	for gene_id in root.findall('gene_ID'):
		gene_name = gene_id.get('value')
		gene_names.add(gene_name)

	gene_names_list = list(gene_names)
	genefinder = pd.DataFrame({'gene_name': gene_names_list, 'source': 'Genefinder'})

	return(genefinder)


# Process BLASTN hits
def extract_blast_hits(blast_results):
	blast_hits = pd.read_csv(blast_results, sep='\t')
	gene_ids = blast_hits['sseqid']
	unique_gene_ids = gene_ids.unique()
	blast = pd.DataFrame({'gene_name': unique_gene_ids, 'source': 'BLASTN'})

	return(blast)


# Process VirulenceFinder hits (optional)
def extract_virulencefinder_hits(virulencefinder_results):
	if virulencefinder_results == "None":
		virulencefinder = pd.DataFrame(columns=['gene_name', 'source'])
	else:
		vf_hits = pd.read_csv(virulencefinder_results, sep='\t')
		gene_ids = vf_hits['Virulence factor']
		unique_gene_ids = gene_ids.unique()
		virulencefinder = pd.DataFrame({'gene_name': unique_gene_ids, 'source': 'VirulenceFinder'})

	return(virulencefinder)


# Merge the results of BLASTN, GeneFinder and (optional) VirulenceFinder
def merge_results(blast, genefinder, virulencefinder):
	original_gene_ids_blast = blast['gene_name'].copy()

	if genefinder.shape[0]>0:
		original_gene_ids_genefinder = genefinder['gene_name'].copy()
		genefinder['gene_name'] = genefinder['gene_name'].str.lower()
	else:
		original_gene_ids_genefinder = pd.Series([], name ='gene_name')
		genefinder['gene_name'] = genefinder['gene_name'].to_string()

	original_gene_ids_vf = virulencefinder['gene_name'].copy()
	original_gene_ids = pd.concat([original_gene_ids_blast, original_gene_ids_genefinder, original_gene_ids_vf]).unique()

	blast['gene_name'] = blast['gene_name'].str.lower()
	virulencefinder['gene_name'] = virulencefinder['gene_name'].str.lower()

	merged_hits = pd.merge(blast, genefinder, on="gene_name", how='outer', suffixes=('_1', '_2'))
	merged_hits = pd.merge(merged_hits, virulencefinder, on="gene_name", how='outer')
	merged_hits['source'] = merged_hits[['source_1', 'source_2', 'source']].apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)

	gene_id_map = {gene_name.lower(): gene_name for gene_name in original_gene_ids}
	merged_hits['gene_name'] = merged_hits['gene_name'].map(gene_id_map).fillna(merged_hits['gene_name'])

	merged_hits = merged_hits[['gene_name', 'source']]

	return merged_hits


# Write the results to long and short summaries
def write_results(merged_hits, bin, output):
	merged_hits.to_csv(output+ '.' + bin +".gene_hits.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')

	collapsed_gene_names = ', '.join(merged_hits['gene_name'])

	hit_summary = pd.DataFrame({'Sample':[output], 'Bin': [bin],'Genes identified': [collapsed_gene_names]})
	hit_summary.to_csv(output+ '.' + bin +".gene_hits.summary.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')


# Parse arguments from the command line.
def parse_args():
	description = 'Merge targeted gene finding results, BLASTN, GeneFinder and VirulenceFinder (optional). Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--blastn', required=True, help="BLASTN results")
	parser.add_argument('--genefinder', required=True, help="GeneFinder results")
	parser.add_argument('--virulencefinder', required=False, help="VirulenceFinder results")
	parser.add_argument('--output', required=True, help="Output TSV name")
	parser.add_argument('--bin', required=True, help="Bin name")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()

	genefinder = extract_genefinder_hits(args.genefinder)
	blast = extract_blast_hits(args.blastn)

	if args.virulencefinder:
		virulencefinder = extract_virulencefinder_hits(args.virulencefinder)
	else:
		virulencefinder = extract_virulencefinder_hits("None")

	merged_hits = merge_results(blast, genefinder, virulencefinder)
	write_results(merged_hits, args.bin, args.output)

if __name__ == "__main__":
        main()
