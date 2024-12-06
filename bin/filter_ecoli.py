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

# Process results from Shigeifinder
def process_shigeifinder(shigeifinder_results):
	shigeifinder_designation = []
	shigeifinder = pd.read_csv(shigeifinder_results, sep='\t')
	shigeifinder.columns = ['Name','Shigeifinder: ipaH', 'Shigeifinder: Virulence_plasmid', 'Shigeifinder: Cluster', 'Shigeifinder: Serotype','Shigeifinder: O-antigen', 'Shigeifinder: H-antigen', 'Shigeifinder: Notes']
	shigeifinder_clusters = ["C1", "C2", "C3", "CSS", "CSB12", "CSB13", "CSB13-atypical", "CSD1", "CSD8", "CSD10"]
	if shigeifinder['Shigeifinder: Cluster'].isin(shigeifinder_clusters).any() == True:
		shigeifinder_designation = "Shigella sp."
	else:
		shigeifinder_designation = ["Unknown"]

	return(shigeifinder_designation, shigeifinder)


# Process results from Shigatyper
def process_shigatyper(shigatyper_results):
	shigatyper_designation = []
	shigatyper = pd.read_csv(shigatyper_results, sep='\t')
	shigatyper.columns = ['Shigatyper: input','Shigatyper: prediction', 'Shigatyper: ipaB', 'Shigatyper: notes']
	shigatyper['Name'] = shigatyper['Shigatyper: input'].str.split('.')[0][0]

	if shigatyper.at[0, 'Shigatyper: prediction'] == "" or shigatyper.at[0, 'Shigatyper: prediction'] == "Not Shigella or EIEC" or shigatyper.at[0, 'Shigatyper: prediction'].startswith("No prediction", "Unable"):
		shigatyper_designation = ["Unknown"]
	elif shigatyper.at[0,'Shigatyper: prediction'].startswith("EIEC"):
		shigatyper_designation = ["EIEC"]
	elif m5.at[0, 'Shigatyper: prediction'].startswith("Shigella boydii"):
		shigatyper_designation = ["Shigella boydii"]
	elif m5.at[0, 'Shigatyper: prediction'].startswith("Shigella dysenteriae"):
		shigatyper_designation = ["Shigella dysenteriae"]
	elif m5.at[0, 'Shigatyper: prediction'].startswith("Shigella flexneri"):
		shigatyper_designation = ["Shigella flexneri"]
	elif m5.at[0, 'Shigatyper: prediction'].startswith("Shigella sonnei"):
		shigatyper_designation = ["Shigella sonnei"]
	else:
		shigatyper_designation = ["Unknown"]

	return(shigatyper_designation, shigatyper)


# Process results from STECfinder (requires filtering based on both cluster and cluster notes)
def process_stecfinder(stecfinder_results):
	stecfinder_designation = []
	stecfinder = pd.read_csv(stecfinder_results, sep='\t')
	stecfinder.columns = ['Name', 'STECFinder: Cluster', 'STECFinder: Cluster Serotype', 'STECFinder: Serotype','STECFinder: Big10 serotype', 'STECFinder: O-antigens','STECFinder: H-antigens', 'STECFinder: stx type', 'STECFinder: ipaH presence', 'STECFinder: Notes']
	stecfinder_clusters = ["AM1","AM10","AM11","AM12","AM13","AM14","AM15","AM16","AM17","AM18","AM19","AM2","AM20","AM21","AM22","AM23","AM24","AM25","AM26","AM27","AM28","AM29","AM3","AM30","AM31","AM32","AM33","AM34","AM35","AM36","AM37","AM4","AM5","AM6","AM7","AM8","AM9","B1M1","B1M10","B1M100","B1M101","B1M102","B1M103","B1M104","B1M105","B1M106","B1M107","B1M108","B1M109","B1M11","B1M110","B1M111","B1M112","B1M113","B1M114","B1M115","B1M116","B1M117","B1M118","B1M119","B1M12","B1M120","B1M121","B1M122","B1M123","B1M124","B1M125","B1M126","B1M13","B1M14","B1M15","B1M16","B1M17","B1M18","B1M19","B1M2","B1M20","B1M21","B1M22","B1M23","B1M24","B1M25","B1M26","B1M27","B1M28","B1M29","B1M3","B1M30","B1M31","B1M32","B1M33","B1M34","B1M35","B1M36","B1M37","B1M38","B1M39","B1M4","B1M40","B1M41","B1M42","B1M43","B1M44","B1M45","B1M46","B1M47","B1M48","B1M49","B1M5","B1M50","B1M51","B1M52","B1M53","B1M54","B1M55","B1M56","B1M57","B1M58","B1M59","B1M6","B1M60","B1M61","B1M62","B1M63","B1M64","B1M65","B1M66","B1M67","B1M68","B1M69","B1M7","B1M70","B1M71","B1M72","B1M73","B1M74","B1M75","B1M76","B1M77","B1M78","B1M79","B1M8","B1M80","B1M81","B1M82","B1M83","B1M84","B1M85","B1M86","B1M87","B1M88","B1M89","B1M9","B1M90","B1M91","B1M92","B1M93","B1M94","B1M95","B1M96","B1M97","B1M98","B1M99","B2M1","B2M10","B2M11","B2M12","B2M13","B2M14","B2M2","B2M3","B2M4","B2M5","B2M6","B2M7","B2M8","B2M9","C1","C10","C11","C12","C13","C14","C15","C16","C17","C18","C2","C3","C4","C5","C6","C7","C8","C9","CM1","CM2","CM3","CM4","CM5","CM6","CM7","DM1","DM10","DM11","DM12","DM13","DM14","DM15","DM16","DM17","DM18","DM19","DM2","DM20","DM21","DM22","DM3","DM4","DM5","DM6","DM7","DM8","DM9","EM1","EM10","EM11","EM12","EM13","EM14","EM15","EM16","EM17","EM18","EM19","EM2","EM3","EM4","EM5","EM6","EM7","EM8","EM9","GM1","GM2","GM3","GM4","O103H11","O103H2","O111H8","O118H16","O123H2","O157H7","O26H11","O45H2","O45H2-AM37","O45H2-C3"]
	stecfinder_notes = ["ipaH+stx+ = Possible EIEC/Shigella, try out other tool shigeifinder!","Strain with cluster specific genes but no stx detected","ipaH+stx- = Possible EIEC/Shigella, try out other tool shigeifinder!","Strain with cluster specific genes but no stx detected and ipaH detected - looks like EIEC/Shigella but clusters with STEC","Strain in STEC cluster but has ipaH and stx - looks like EIEC/Shigella but clusters with STEC", "Unexpected genetic combination, possible bug, please post issue on shigeifinder github repo", " Missmatch between big10 specific genes and antigen genes, possible big10 specific genes false positive.", " Possible missing antigen in big10 isolate."," Detected O and H antigen genes do not match previously known serotypes for this cluster"]

	if stecfinder['STECFinder: Cluster'].isin(stecfinder_clusters).any() == True:
		for val in stecfinder_notes:
			if stecfinder['STECFinder: Notes'].astype(str).str.contains(val).any() == True:
				stecfinder_designation = ["Unknown"]
				break
			elif stecfinder.at[0, 'STECFinder: Notes'] == "ipaH-stx- = Non-STEC E.coli":
				stecfinder_designation = ["Non-STEC E.coli"]
			elif stecfinder.at[0, 'STECFinder: Notes'] == "STEC not from any major STEC lineages.":
				stecfinder_designation = ["STEC"]
			else:
				stecfinder_designation = ["STEC"]

	elif stecfinder.at[0, 'STECFinder: Cluster'] == "Unclustered STEC":
		stecfinder_designation = ["STEC"]
	elif stecfinder.at[0, 'STECFinder: Cluster'].startswith("EIEC/Shigella") or stecfinder.at[0, 'STECFinder: Cluster'] == "-":
		stecfinder_designation = ["Unknown"]
	elif stecfinder.at[0, 'STECFinder: Cluster'].startswith("Other_Ecoli"):
		stecfinder_designation = ["Non-STEC E.coli"]
	else:
		stecfinder_designation = ["Unknown"]

	return(stecfinder_designation, stecfinder)

# Process results from ECtyper
def process_ectyper(ectyper_results):
	ectyper_designation = []
	ectyper = pd.read_csv(ectyper_results, sep='\t')
	ectyper.columns = ['Ectyper: input','Ectyper: Species','Ectyper: O-type','Ectyper: H-type','Ectyper: Serotype','Ectyper: QC','Ectyper: Evidence','Ectyper: GeneScores','Ectyper: AlleleKeys','Ectyper: GeneIdentities(%)','Ectyper: GeneCoverages(%)','Ectyper: GeneContigNames','Ectyper: GeneRanges','Ectyper: GeneLengths','Ectyper: Database','Ectyper: Warnings']
	ectyper['Name'] = ectyper['Ectyper: input'].str.split('.')[0][0]

	if ectyper.at[0, 'Ectyper: Species'].startswith("Shigella flexneri"):
		ectyper_designation = ["Shigella flexneri"]
	elif ectyper.at[0, 'Ectyper: Species'].startswith("Shigella sonnei"):
		ectyper_designation = ["Shigella sonnei"]
	elif ectyper.at[0, 'Ectyper: Species'].startswith("Shigella dysenteriae"):
		ectyper_designation = ["Shigella dysenteriae"]
	elif ectyper.at[0, 'Ectyper: Species'].startswith("Shigella boydii"):
		ectyper_designation = ["Shigella boydii"]
	elif ectyper.at[0, 'Ectyper: Species'].startswith("Escherichia coli"):
		ectyper_designation = ["Non-STEC E.coli"]
	else:
		ectyper_designation = ["Unknown"]

	return(ectyper_designation, ectyper)

# Process results from mlst and Krocus (id'ing species based on known associations between sequence types of pathotypes)
def process_mlst(mlst_results, krocus_results):
	mlst_designation = []
	mlst = pd.read_csv(mlst_results, sep='\t')
	mlst = mlst[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
	mlst.columns = ['Name', "MLST: Sequence type", "Internal UKHSA: Clonal complex (mlst)"]
	krocus = pd.read_csv(krocus_results, sep='\t')
	krocus = krocus[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
	krocus.columns = ['Name', "Krocus: Sequence type", "Internal UKHSA: Clonal complex (Krocus)"]

	cc_merged = pd.merge(mlst, krocus,on=['Name'], how = 'outer')
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

	if cc_merged.at[0, 'Internal UKHSA: Clonal complex'] == "CC245":
		mlst_designation = ["Shigella flexneri"]
	elif cc_merged.at[0, 'Internal UKHSA: Clonal complex'] == "CC152":
		mlst_designation = ["Shigella sonnei"]
	elif cc_merged.at[0, 'Internal UKHSA: Clonal complex'] == "CC145" or cc_merged.at[0, 'Internal UKHSA: Clonal complex'] == "CC288":
		mlst_designation = ["Shigella sp."]
	elif cc_merged.at[0, 'Sequence type'] == 148:
		mlst_designation = ["Shigella dysenteriae"]
	else:
		mlst_designation = ["Unknown"]

	return(mlst_designation, cc_merged)


# Process results from Mykrobe (only relevant for Shigella sonnei)
def process_mykrobe(mykrobe_results):
	mykrobe_designation = []
	if mykrobe_results == "None":
		mykrobe_designation = ["Unknown"]
		mykrobe = pd.DataFrame(columns=['Name','Mykrobe: input', 'Mykrobe: genotype_model', 'Mykrobe: kmer_size', 'Mykrobe: phylo_group', 'Mykrobe: species', 'Mykrobe: lineage', 'Mykrobe: phylo_group_per_covg','Mykrobe: species_per_covg', 'Mykrobe: lineage_per_covg','Mykrobe: phylo_group_depth', 'Mykrobe: species_depth', 'Mykrobe: lineage_depth'], index=range(1))
	else:
		mykrobe = pd.read_csv(mykrobe_results)
		mykrobe = mykrobe[['sample','genotype_model','kmer_size','phylo_group','species','lineage','phylo_group_per_covg','species_per_covg','lineage_per_covg','phylo_group_depth','species_depth','lineage_depth']]
		mykrobe.columns = ['Mykrobe: input', 'Mykrobe: genotype_model', 'Mykrobe: kmer_size', 'Mykrobe: phylo_group', 'Mykrobe: species', 'Mykrobe: lineage', 'Mykrobe: phylo_group_per_covg','Mykrobe: species_per_covg', 'Mykrobe: lineage_per_covg','Mykrobe: phylo_group_depth', 'Mykrobe: species_depth', 'Mykrobe: lineage_depth']
		mykrobe['Name'] = mykrobe['Mykrobe: input'].str.split('.')[0][0]

		if mykrobe.at[0, 'Mykrobe: species'] == "Shigella_sonnei":
			mykrobe_designation = ["Shigella sonnei"]
		else:
			mykrobe_designation = ["Unknown"]

	return(mykrobe_designation, mykrobe)


# Create a consensus pathotype by combining all available evidence
def consensus_designation(shigeifinder_designation, shigatyper_designation, stecfinder_designation, mlst_designation, ectyper_designation, mykrobe_designation):
	merged_degs = pd.DataFrame({'tag':shigeifinder_designation + shigatyper_designation + stecfinder_designation + mlst_designation +ectyper_designation + mykrobe_designation})
	merged_degs = merged_degs.value_counts().reset_index(name='count')
	merged_degs = merged_degs[(merged_degs["count"]> 1) & (merged_degs["tag"] != "Unknown")]
	if merged_degs.empty:
		top_assignment = "Unknown"

	else:
		merged_degs = merged_degs.groupby('count').agg('; '.join).reset_index()
		merged_degs = merged_degs.sort_values(by='count', ascending=False).iloc[0]
		id_string = merged_degs['tag'] + ' (' + merged_degs['count'].astype(str) + ')'
		top_assignment = (id_string).replace('_count','').replace('_',' ')

	return(top_assignment)


# Merge metrics into a single tables returning the relevant columns
def merge_tables(shigeifinder, shigatyper, stecfinder, mlst, ectyper, mykrobe, top_assignment, output, bin):
	metrics = [shigeifinder, shigatyper, stecfinder, ectyper, mlst, mykrobe]
	mykrobe['Name'] = shigeifinder['Name'][0]        
	shigatyper['Name'] = shigeifinder['Name'][0]
	merge = partial(pd.merge, on=['Name'], how='outer')
	merged_metrics = (reduce(merge, metrics))

	merged_metrics['Inferred classification'] = top_assignment
	merged_metrics['Name'] = output
	merged_metrics['Bin'] = bin
	merged_metrics = merged_metrics[['Name','Bin','Inferred classification','Shigeifinder: ipaH','Shigeifinder: Virulence_plasmid',
			 'Shigeifinder: Cluster','Shigeifinder: Serotype','Shigeifinder: O-antigen',
			 'Shigeifinder: H-antigen','Shigeifinder: Notes','Shigatyper: prediction',
			 'Shigatyper: ipaB','Shigatyper: notes','STECFinder: Cluster','STECFinder: Cluster Serotype',
			 'STECFinder: Serotype','STECFinder: Big10 serotype','STECFinder: O-antigens','STECFinder: H-antigens',
			 'STECFinder: stx type','STECFinder: ipaH presence','STECFinder: Notes','Ectyper: Species',
			 'Ectyper: O-type','Ectyper: H-type','Ectyper: Serotype','Ectyper: Evidence',
			 'Ectyper: AlleleKeys','Sequence type','Internal UKHSA: Clonal complex','Mykrobe: phylo_group',
			 'Mykrobe: species','Mykrobe: lineage']]

	merged_metrics.to_csv(output+".ecoli_typing_detailed.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')

	return(merged_metrics)


# Simplify merged table
def simplify_tables(merged_metrics, output):
	shigeifinder_serotype = merged_metrics.agg('{0[Shigeifinder: O-antigen]};{0[Shigeifinder: H-antigen]}'.format, axis=1)[0]
	if shigeifinder_serotype != "nan;nan":
		xlist = [shigeifinder_serotype, merged_metrics['STECFinder: Cluster Serotype'][0], merged_metrics['STECFinder: Serotype'][0],merged_metrics['STECFinder: Big10 serotype'][0], merged_metrics['Ectyper: Serotype'][0]]
	else:
		xlist = [merged_metrics['STECFinder: Cluster Serotype'][0], merged_metrics['STECFinder: Serotype'][0],merged_metrics['STECFinder: Big10 serotype'][0], merged_metrics['Ectyper: Serotype'][0]]

	ylist = [merged_metrics['Shigeifinder: ipaH'][0],merged_metrics['STECFinder: ipaH presence'][0]]
	x1list = pd.Series(xlist).dropna().drop_duplicates().tolist()

	x1list = (', '.join(map(str, (set(x1list)))))
	y1list = (', '.join(map(str, sorted(set(ylist)))))
	merged_metrics['Consensus serotype'] = x1list
	merged_metrics['Consensus ipaH'] = y1list

	merged_metrics_short = (merged_metrics[['Name','Bin','Inferred classification','Ectyper: Species','Mykrobe: species','Mykrobe: lineage','Mykrobe: phylo_group',
						'Shigeifinder: Cluster','STECFinder: Cluster','Consensus serotype', 'Consensus ipaH',
						'Shigatyper: ipaB','Shigeifinder: Virulence_plasmid','STECFinder: stx type',
						'Sequence type','Internal UKHSA: Clonal complex','Shigatyper: prediction']])

	merged_metrics_short.to_csv(output+".ecoli_typing_short.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')


# Parse arguments from the command line.
def parse_args():
	description = 'Identify most likely species/pathotype identity for input metagenome assembled genome. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--shigeifinder', required=True, help="Shigeifinder results")
	parser.add_argument('--shigatyper', required=True, help="Shigatyper results")
	parser.add_argument('--ectyper', required=True, help="ECtyper results")
	parser.add_argument('--stecfinder', required=True, help="STECFinder results")
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
	shigeifinder_designation, shigeifinder = process_shigeifinder(args.shigeifinder)
	shigatyper_designation, shigatyper = process_shigatyper(args.shigatyper)
	stecfinder_designation, stecfinder = process_stecfinder(args.stecfinder)
	mlst_designation, mlst = process_mlst(args.mlst, args.krocus)
	ectyper_designation, ectyper = process_ectyper(args.ectyper)
	if args.mykrobe:
		mykrobe_designation, mykrobe = process_mykrobe(args.mykrobe)
	else:
		mykrobe_designation, mykrobe = process_mykrobe("None")
	top_assignment = consensus_designation(shigeifinder_designation, shigatyper_designation, stecfinder_designation, mlst_designation, ectyper_designation, mykrobe_designation)
	merged_metrics = merge_tables(shigeifinder, shigatyper, stecfinder, mlst, ectyper, mykrobe, top_assignment, args.output, args.bin)
	simplify_tables(merged_metrics, args.output)

if __name__ == "__main__":
	main()
