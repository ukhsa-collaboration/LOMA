#!/usr/bin/env python

# Filter and merge Kraken2, Centrifuger or Sylph results against a database of target species

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports
import pandas as pd
import numpy as np
import argparse


###### Functions

# Process Kraken2 or Centrifuger hits
def process_hits(taxhits, target_list, min_reads, min_frac):
	taxhits_tbl = pd.read_csv(taxhits, delimiter='\t')
	targets_tbl = pd.read_csv(target_list, delimiter='\t', header=None)
	targets_tbl.columns = ['Species', 'taxonomy_id', 'Kingdom', 'Rank']
	targets_tbl.dropna(subset=["taxonomy_id"], inplace=True)
	targets_tbl = targets_tbl.drop_duplicates(subset='taxonomy_id', keep='first')

	merged_taxhits = pd.merge(taxhits_tbl, targets_tbl, on=['taxonomy_id'])

	merged_taxhits_x = merged_taxhits[['Species', 'Rank', 'Kingdom', 'taxonomy_id', 'name','new_est_reads', 'fraction_total_reads']]

	merged_taxhits_x = merged_taxhits_x[(merged_taxhits_x.new_est_reads > int(min_reads)) & (merged_taxhits_x.fraction_total_reads > min_frac)]
	merged_union = merged_taxhits_x
	merged_union.columns = ['Species', 'Rank', 'Kingdom','Taxonomy ID', 'Name', 'Estimated reads', 'Fraction total reads' ]
	merged_union = merged_union.drop('Name', axis=1)

	merged_union.fillna('', inplace=True)
	merged_union['Estimated reads'] = merged_union['Estimated reads'].astype('Int64')

	return merged_union


def process_sylph(taxhits, target_list, gtdbtk_fn, min_frac):
	syl_fn_df = pd.read_csv(gtdbtk_fn, compression='gzip', sep='\t', header=None)
	syl_fn_df.columns = ['Genome_file_1','tax']

	targets_tbl = pd.read_csv(target_list, delimiter='\t', header=None)
	targets_tbl.columns = ['Species', 'taxonomy_id', 'Kingdom', 'Rank']

	syl = pd.read_csv(taxhits, delimiter='\t')
	syl['Genome_file_1'] = syl['Genome_file'].str.split('/').str[-1].str.split("_").str[0:2].apply('_'.join)

	syl = syl[(syl['Taxonomic_abundance'] > (min_frac*100))]

	syl_merge = pd.merge(syl,syl_fn_df,on=['Genome_file_1'], how = 'left')
	syl_merge['Species1'] = syl_merge['tax'].str.replace(r'd__;p__;o__;f__;g__;s__$', '', regex=True).str.replace(r';p__;o__;f__;g__;s__$', '', regex=True).str.replace(r';c__;o__;f__;g__;s__$', '', regex=True).str.replace(r';o__;f__;g__;s__$', '', regex=True).str.replace(r';f__;g__;s__$', '', regex=True).str.replace(r';g__;s__$', '', regex=True).str.replace(r';s__$', '', regex=True).str.split(';').str[-1]

	syl_merge['rank'] = np.select(
		[syl_merge['Species1'].str.contains('s__'),
		syl_merge['Species1'].str.contains('g__'),
		syl_merge['Species1'].str.contains('f__'),
		syl_merge['Species1'].str.contains('o__'),
		syl_merge['Species1'].str.contains('c__'),
		syl_merge['Species1'].str.contains('p__'),
		syl_merge['Species1'].str.contains('d__'),
		syl_merge['Species1'].str.contains('Unclassified')],
		['Species', 'Genus','Family','Order','Class','Phylum','Domain','Unclassified'],
		default=np.nan,
	)
	syl_merge['Species'] = syl_merge['Species1'].str.replace('s__', '').str.replace('g__', '').str.replace('f__', '').str.replace('o__', '').str.replace('c__', '').str.replace('p__', '').str.replace('d__', '')
	syl_merge = syl_merge[["Taxonomic_abundance", "Sequence_abundance","Species"]]

	syl_merged_taxhits = pd.merge(syl_merge, targets_tbl, on=['Species'])

	syl_merged_taxhits = syl_merged_taxhits[["Species", 'Rank', "Kingdom", "taxonomy_id", "Taxonomic_abundance", "Sequence_abundance"]]
	syl_merged_taxhits.columns = ["Species", 'Rank', "Kingdom", "Taxonomy ID", "Taxonomic abundance", "Sequence abundance"]

	syl_merged_taxhits['Taxonomy ID'] = syl_merged_taxhits['Taxonomy ID'].astype('Int64')

	return(syl_merged_taxhits)

# Save results to output file
def write_hits_to_file(taxhits, merged_union, output):
	tool = taxhits.split(".")[-2]
	output_file = output + '.' + tool + '.target_species.tsv'
	merged_union.to_csv(output_file, sep='\t', index=False)


# Parse arguments from the command line
def parse_args():
	description = 'Filter Kraken2, Centrifuger or Sylph results against a database of target species. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--taxhits', required=True, help="Results from profiling tool (Either Kraken2 or Centrifuger)")
	parser.add_argument('--output', required=True, help="Output TSV name")
	parser.add_argument('--targets', required=True, help="List of target species")
	parser.add_argument('--min_reads', required=True, help="Minimum read count to include hit")
	parser.add_argument('--min_frac', required=True, type=float, help="Minimum read fraction to include hit")
	parser.add_argument('--mode', required=True, help="Determines how input is processed [Bracken|Sylph]")
	parser.add_argument('--gtdb_fn', type=str, required=False, help='Path to gtdb_*_metadata.tsv.gz file')
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()

	if args.mode == "Bracken":
		merged_union = process_hits(args.taxhits, args.targets, args.min_reads, args.min_frac)
		write_hits_to_file(args.taxhits, merged_union, args.output)

	elif args.mode == "Sylph" and args.gtdb_fn:
		syl_merged = process_sylph(args.taxhits, args.targets, args.gtdb_fn, args.min_frac)
		write_hits_to_file(args.taxhits, syl_merged, args.output)

	elif args.mode == "Sylph" and not args.gtdb_fn:
		print("--mode Sylph specific but gtdb_*_metadata.tsv.gz was not provided, please specify using the --gtdb_fn parameter")
		exit
	else:
		print("Incorrect input for parameter --mode , it should be either 'Bracken' or 'Sylph'")
		exit

if __name__ == "__main__":
	main()
