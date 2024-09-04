#!/usr/bin/env python

# Merges results from multiple E. coli/Shigella typing tools to provide a summarized output and a estimation of assignment/species/pathotype.
# Note: Hard-coded cluter IDs in places due to the way various tools function. These may need to be updated at some point.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports
import pandas as pd
import numpy as np
import argparse
from functools import partial, reduce

###### Functions

# Process results from SeqSero2
def process_seqsero2(seqsero2_results, bin):
	seqsero2 = pd.read_csv(seqsero2_results, sep='\t')
	seqsero2['Name'] = seqsero2['Sample name'].str.split('.')[0][0]
	seqsero2['Bin'] = bin
	seqsero2 = seqsero2[["Name", "Bin", "Predicted antigenic profile", "Predicted identification", "Predicted serotype", "Note"]]
	seqsero2.columns = ["Name", "Bin", "SeqSero2: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))", "Seqsero2: Predicted identification","SeqSero2: Predicted serotype", "SeqSero2: Note"]

	return(seqsero2)


# Process results from SISTR
def process_sistr(sistr_results, bin):
	sistr = pd.read_csv(sistr_results, sep='\t')
	sistr['Name'] = sistr['genome'].str.split('.')[0][0]
	sistr['Bin'] = bin
	sistr["sero_temp"] = str(sistr["o_antigen"][0] + ":" + sistr["h1"][0] + ":" + sistr["h2"][0])
	sistr = sistr[["Name", "Bin", "cgmlst_subspecies", "cgmlst_ST","qc_messages", "sero_temp","serogroup","serovar", "serovar_antigen","serovar_cgmlst" ]]
	sistr.columns = ["Name", "Bin", "SISTR: cgMLST subspecies", "SISTR: cgMLST sequence type", "SISTR: QC messages", "SISTR: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))","SISTR: Serogroup", "SISTR: Serovar", "SISTR: Serovar antigen", "SISTR: Serovar cgMLST"]

	return(sistr)


# Process results from mlst and Krocus (id'ing species based on known associations between sequence types of pathotypes)
def process_mlst(mlst_results, krocus_results, bin):
	mlst = pd.read_csv(mlst_results, sep='\t')
	mlst = mlst[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
	mlst.columns = ['Name', "MLST: Sequence type", "Internal UKHSA: Clonal complex (mlst)"]
	krocus = pd.read_csv(krocus_results, sep='\t')
	krocus = krocus[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
	krocus.columns = ['Name', "Krocus: Sequence type", "Internal UKHSA: Clonal complex (Krocus)"]

	cc_merged = pd.merge(mlst, krocus,on=['Name'], how = 'outer')
	cc_merged['Bin'] = bin
	cc_merged = cc_merged.fillna("NA")

	if cc_merged.at[0, 'Internal UKHSA: Clonal complex (mlst)'] == cc_merged.at[0, 'Internal UKHSA: Clonal complex (Krocus)']:
		cc_merged['Internal UKHSA: Clonal complex'] = cc_merged.at[0, 'Internal UKHSA: Clonal complex (mlst)']
	elif cc_merged.at[0, 'Internal UKHSA: Clonal complex (mlst)'] == "NA" and cc_merged.at[0, 'Internal UKHSA: Clonal complex (Krocus)'] != "NA":
		cc_merged['Internal UKHSA: Clonal complex'] = cc_merged.at[0, 'Internal UKHSA: Clonal complex (Krocus)']
	elif cc_merged.at[0, 'Internal UKHSA: Clonal complex (mlst)'] != "NA" and cc_merged.at[0, 'Internal UKHSA: Clonal complex (Krocus)'] == "NA":
		cc_merged['Internal UKHSA: Clonal complex'] = cc_merged.at[0, 'Internal UKHSA: Clonal complex (mlst)']
	else:
		cc_merged['Internal UKHSA: Clonal complex'] = "Unknown"

	if cc_merged.at[0, 'MLST: Sequence type'] == cc_merged.at[0, 'Krocus: Sequence type']:
		cc_merged['Sequence type'] = cc_merged.at[0, 'MLST: Sequence type']
	elif cc_merged.at[0, 'MLST: Sequence type'] == "-" and cc_merged.at[0, 'Krocus: Sequence type'] != "ND":
		cc_merged['Sequence type'] = cc_merged.at[0, 'Krocus: Sequence type']
	elif cc_merged.at[0, 'MLST: Sequence type'] != "-" and cc_merged.at[0, 'Krocus: Sequence type'] == "ND":
		cc_merged['Sequence type'] = cc_merged.at[0, 'MLST: Sequence type']
	else:
		cc_merged['Sequence type'] = "Unknown"

	return(cc_merged)


# Process results from Mykrobe
def process_mykrobe(mykrobe_results, bin):
	if mykrobe_results == "None":
		mykrobe = pd.DataFrame(columns=['Name','Bin', 'Mykrobe: input', 'Mykrobe: genotype_model', 'Mykrobe: kmer_size', 'Mykrobe: phylo_group', 'Mykrobe: species', 'Mykrobe: lineage', 'Mykrobe: phylo_group_per_covg','Mykrobe: species_per_covg', 'Mykrobe: lineage_per_covg','Mykrobe: phylo_group_depth', 'Mykrobe: species_depth', 'Mykrobe: lineage_depth'], index=range(1))
	else:
		mykrobe = pd.read_csv(mykrobe_results)
		mykrobe = mykrobe[['sample','genotype_model','kmer_size','phylo_group','species','lineage','phylo_group_per_covg','species_per_covg','lineage_per_covg','phylo_group_depth','species_depth','lineage_depth']]
		mykrobe.columns = ['Mykrobe: input', 'Mykrobe: genotype_model', 'Mykrobe: kmer_size', 'Mykrobe: phylo_group', 'Mykrobe: species', 'Mykrobe: lineage', 'Mykrobe: phylo_group_per_covg','Mykrobe: species_per_covg', 'Mykrobe: lineage_per_covg','Mykrobe: phylo_group_depth', 'Mykrobe: species_depth', 'Mykrobe: lineage_depth']
#		mykrobe['Name'] = mykrobe['Mykrobe: input'].str.split('.')[0][0]

	return(mykrobe)


# Merge metrics into a single tables returning the relevant columns
def merge_tables(seqsero2, sistr, mlst, mykrobe, output, bin):
	metrics = [seqsero2, sistr, mlst, mykrobe]
	mykrobe["Name"] = seqsero2["Name"]
	mykrobe["Bin"] = bin
	merge = partial(pd.merge, on=['Name', 'Bin'], how='outer')
	merged_metrics = (reduce(merge, metrics))

	merged_metrics = merged_metrics[['Name','Bin',
					"Seqsero2: Predicted identification", "SeqSero2: Predicted serotype", "SeqSero2: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))",
					"SeqSero2: Note","SISTR: cgMLST subspecies",
					"SISTR: cgMLST sequence type","SISTR: QC messages","SISTR: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))",
					"SISTR: Serogroup","SISTR: Serovar","SISTR: Serovar antigen","SISTR: Serovar cgMLST",'Sequence type','Internal UKHSA: Clonal complex',
					'Mykrobe: phylo_group','Mykrobe: species','Mykrobe: lineage']]

	merged_metrics.to_csv(output+".salmonella_typing_detailed.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')

	return(merged_metrics)


# Simplify merged table
def simplify_tables(merged_metrics, output):
	merged_metrics_short = (merged_metrics[['Name','Bin', "Seqsero2: Predicted identification", "SeqSero2: Predicted serotype",
						"SeqSero2: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))",
						"SISTR: cgMLST subspecies",
						"SISTR: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))",
						"SISTR: Serogroup","SISTR: Serovar","SISTR: Serovar antigen","SISTR: Serovar cgMLST",
						'Sequence type','Internal UKHSA: Clonal complex',
						'Mykrobe: species','Mykrobe: lineage','Mykrobe: phylo_group']])

	merged_metrics_short.to_csv(output+".salmonella_typing_short.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')


# Parse arguments from the command line.
def parse_args():
	description = 'Identify most likely species/pathotype identity for input metagenome assembled genome. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--seqsero2', required=True, help="SeqSero2 results")
	parser.add_argument('--sistr', required=True, help="SISTR results")
	parser.add_argument('--mykrobe', required=False, help="Mykrobe results")
	parser.add_argument('--mlst', required=True, help="mlst results")
	parser.add_argument('--krocus', required=True, help="Krocus results")
	parser.add_argument('--output', required=True, help="Output TSV name")
	parser.add_argument('--bin', required=True, help="Bin name")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()
	mlst = process_mlst(args.mlst, args.krocus, args.bin)
	seqsero2 = process_seqsero2(args.seqsero2, args.bin)
	sistr = process_sistr(args.sistr, args.bin)
	if args.mykrobe:
		mykrobe = process_mykrobe(args.mykrobe, args.bin)
	else:
		mykrobe = process_mykrobe("None", args.bin)
	merged_metrics = merge_tables(seqsero2, sistr, mlst, mykrobe, args.output, args.bin)
	simplify_tables(merged_metrics, args.output)

if __name__ == "__main__":
	main()
