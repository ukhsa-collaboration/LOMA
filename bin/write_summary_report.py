#!/usr/bin/env python

# Create a summary HTML report based on the available outputs from a compeleted pipeline run.

__version__ = '0.1'
__date__ = '03-06-2024'
__author__ = 'D.J.BERGER'


###### Imports
import pandas as pd
import numpy as np
import argparse
import datetime
from jinja2 import Environment, FileSystemLoader
import os
import base64

###### Functions


# Create output for pipeline metadata section
def process_pipeline_metadata(pipeline_version):
	current_time = datetime.datetime.now().strftime("%x ; %X")
	run_meta_data = {
		'date_time' : current_time,
		'pipeline_version' : pipeline_version,
	}

	return run_meta_data


# Combine sample metadata
def process_metadata(sample_id, run_id, barcode, sample_type, logo):
	with open(logo, "rb") as image_file:
		sample_data = {
			'sample_id' : sample_id,
			'run_id' : run_id,
			'barcode' : barcode,
			'sample_type' : sample_type,
			'logo_url' : base64.b64encode(image_file.read()).decode('utf-8')
		}

	return sample_data


# Combine read QC statistics
def process_readqc(read_qc_stats, hostreads_list):
	readqc = pd.read_csv(read_qc_stats, sep='\t')
	hostreads = pd.read_csv(hostreads_list, sep='\t', header=None)
	readqc_data = {
		'total_reads_pre_qc' : int(readqc['Pre-QC'].iloc[0]),
		'total_reads_post_qc' : int(readqc['Post-QC'].iloc[0]),
		'total_human_reads' : int(len(hostreads)),
		'pre_read_n50_bp' : int(readqc['Pre-QC'].iloc[4]),
		'read_n50_bp' : int(readqc['Post-QC'].iloc[4]),
		'total_bases_mb' : round((int(readqc['Post-QC'].iloc[7])),2) ,
		'total_bases_pre_qc' : round((int(readqc['Pre-QC'].iloc[7])),2) ,
		'median_read_length_bp' : int(readqc['Post-QC'].iloc[2]),
		'median_read_quality_x' : '{:.2f}'.format(readqc['Post-QC'].iloc[6]),
	}

	return readqc_data


# Process the Bracken (targets) file
def process_bracken_targets(bracken_targets_ccvals, min_read_count, min_read_prop):
	bracken_targets = pd.read_csv(bracken_targets_ccvals, sep='\t')
	bracken_targets['Perc total reads'] = round((bracken_targets['Fraction total reads']*100),2)
	bracken_targets = bracken_targets[bracken_targets['Perc total reads'] >= float(min_read_prop)]
	bracken_targets = bracken_targets[bracken_targets['Estimated reads'] >= int(min_read_count)]
	bracken_targets = bracken_targets[['Species', 'Rank', 'Estimated reads', 'Perc total reads']]

	bracken_targets.columns = ['species', 'Rank', 'counts', 'percentage']

	bracken_targets = bracken_targets.sort_values(by = 'counts', ascending=False)
	bracken_targets_data = bracken_targets.to_dict(orient='records')

	return bracken_targets_data


# Process the Bracken (all hits) file
def process_bracken_all(bracken_all_results, min_read_count, min_read_prop):
	bracken_all = pd.read_csv(bracken_all_results, sep='\t')
	bracken_all['Perc total reads'] = round((bracken_all['fraction_total_reads']*100),2)
	bracken_all = bracken_all[bracken_all['Perc total reads'] >= float(min_read_prop)]
	bracken_all = bracken_all[bracken_all['new_est_reads'] >= int(min_read_count)]
	bracken_all = bracken_all[['name','taxonomy_lvl', 'new_est_reads', 'Perc total reads']]
	bracken_all['taxonomy_lvl'] = bracken_all['taxonomy_lvl'].replace('S','Species')

	bracken_all.columns = ['species', 'tax_lvl', 'counts', 'percentage']

	bracken_all = bracken_all.sort_values(by = 'counts', ascending=False)

	bracken_all_data = bracken_all.to_dict(orient='records')

	return bracken_all_data


# Combine contig summary states and calculate contig N50
def process_contig_summary(contig_summary_ccvals):
	contig_summary = pd.read_csv(contig_summary_ccvals, sep='\t')
	total_length = contig_summary['Contig length (bp)'].sum()
	count_total_contigs = len(contig_summary)

	ctax = contig_summary.dropna(subset=['skani: Contig average nucleotide identity to reference'])
	contig_tax = round(((len(ctax[ctax['skani: Inferred taxononomic classification rank'].str.match('Species')])/count_total_contigs)*100),2)

	filter = contig_summary['Contig name'].str.contains('unbinned')
	length_binned = contig_summary[~filter]['Contig length (bp)'].sum()
	count_binned = round(((len(contig_summary[~filter]['Contig length (bp)'])/count_total_contigs)*100),2)
	perc_binned = round(((length_binned/total_length)*100),2)

	contig_lengths = sorted((contig_summary['Contig length (bp)'].astype(int).tolist()), reverse=True)

	cumulative_length = 0
	half_len = (total_length/2)
	for length in contig_lengths:
		cumulative_length += length
		if cumulative_length >= half_len:
			c50 = length
			break

	contig_data = {
		'total_length' : round((total_length/1000000),2),
		'perc_binned' : perc_binned,
		'count_total_contigs' : count_total_contigs,
		'contig_n50' : c50,
		'count_binned' : count_binned,
		'count_with_tax' : contig_tax
	}

	return contig_data


# Summarize binning statistics
def process_bin_summary(bin_summary_ccvals):
	bin_summary = pd.read_csv(bin_summary_ccvals, sep='\t')
	bin_summary['Bin'] = bin_summary['Bin'].str.replace('bin_', '')
	bin_count = len(bin_summary)

	btaxin = bin_summary.dropna(subset=['Inferred classification rank'])
	xin = len(bin_summary[bin_summary['Inferred classification rank'].isna()])
	btax_spcount = round(((len(btaxin[btaxin['Inferred classification rank'].str.startswith('Species')])/bin_count)*100),2)
	btax_gencount = round(((len(btaxin[btaxin['Inferred classification rank'].str.startswith('Genus')])/bin_count)*100),2)
	btax_famcount = round(((len(btaxin[btaxin['Inferred classification rank'].str.startswith('Family')])/bin_count)*100),2)
	btax_unccount = round((((len(btaxin[btaxin['Inferred classification rank'].str.startswith('Unclassified')])+xin)/bin_count)*100),2)

	bin_fail = round(((len(bin_summary[bin_summary['Bin quality classification'] == "QC fail"])/bin_count)*100),2)
	bin_mq = round(((len(bin_summary[bin_summary['Bin quality classification'] == "Medium quality"])/bin_count)*100),2)
	bin_pa = round(((len(bin_summary[bin_summary['Bin quality classification'] == "Partial assembly"])/bin_count)*100),2)
	bin_hq = round(((len(bin_summary[bin_summary['Bin quality classification'] == "High quality"])/bin_count)*100),2)


	bin_data = {
		'btax_sp' : btax_spcount,
		'btax_gen' : btax_gencount,
		'btax_fam' : btax_famcount,
		'btax_unc' : btax_unccount,
		'bin_count' : bin_count,
		'bin_qc_fail' : bin_fail,
		'bin_qc_hq' : bin_hq,
		'bin_qc_mq' : bin_mq,
		'bin_qc_pa' : bin_pa,
		'bin_qc_fail' : bin_fail
	}

	return(bin_data, bin_summary)


# Subset per-bin statistics
def process_bin_typing(bin_summary):
	bin_summary['Assembly sizes (Mb)'] = round((bin_summary['Assembly length (bp)']/1000000),2)

	bin_summary['Contig count'] = bin_summary['Contigs (count)'].astype(int)
	bin_summary = bin_summary[['Bin', 'Inferred taxononomic classification', 'Inferred classification rank', 'GC (%)','Assembly sizes (Mb)', 'Contig count', 'CheckM: Completeness','CheckM: Contamination','Bin quality classification']]

	bin_summary['Inferred taxononomic classification'] = bin_summary['Inferred taxononomic classification'].fillna("Unclassified")
	bin_summary['Inferred classification rank'] = bin_summary['Inferred classification rank'].fillna("Unclassified")

	bin_summary.columns = ['Bin', 'itc', 'icr', 'gc', 'as', 'cc','checkmCM', 'checkmCT', 'bq']
	bin_summary = bin_summary.sort_values('Bin')
	bin_summary_clean = bin_summary
	bin_summary_data = bin_summary.to_dict(orient='records')

	return(bin_summary_data, bin_summary_clean)


# Process MLST results
def process_mlst_indv(mlst_filenames):
	mlst = pd.concat([pd.read_csv(f, sep='\t') for f in mlst_filenames ])
	mlst = mlst[['Name','Sample bin', "mlst: Scheme", "Sequence type", "Internal UKHSA: Clonal complex"]]
	mlst.columns = ['Name', 'Bin', "scheme", "st", "cc"]

	mlst.fillna("Unknown", inplace=True)
	mlst['st'] = mlst['st'].apply(lambda x: str(x).replace('-','Unknown'))
	mlst['st'] = mlst['st'].str.replace('.0', '')
	mlst['Bin'] = mlst['Bin'].str.replace('bin_', '')

	mlst = mlst.sort_values('Bin')
	cc_merged_data = mlst.to_dict(orient='records')

	return cc_merged_data


def process_krocus_indv(krocus_filenames):
	krocus = pd.concat([pd.read_csv(f, sep='\t') for f in krocus_filenames ])
	krocus = krocus[['Name','Sample bin',"krocus: Scheme", "Sequence type", "Internal UKHSA: Clonal complex"]]
	krocus.columns = ['Name', 'Bin', "scheme", "st", "cc"]

	krocus.fillna("Unknown", inplace=True)
	krocus['st'] = krocus['st'].replace('ND', 'Unknown').replace('.0', '')
	krocus['Bin'] = krocus['Bin'].str.replace('bin_', '')
	krocus = krocus.sort_values('Bin')

	cc_merged_data = krocus.to_dict(orient='records')

	return cc_merged_data



# Merge MLST and Krocus results
def process_mlst(krocus_filenames, mlst_filenames):
	krocus = pd.concat([pd.read_csv(f, sep='\t') for f in krocus_filenames ])
	mlst = pd.concat([pd.read_csv(f, sep='\t') for f in mlst_filenames ])

	mlst = mlst[['Name','Sample bin', "mlst: Scheme", "Sequence type", "Internal UKHSA: Clonal complex"]]
	mlst.columns = ['Name', 'Bin', "mlst: Scheme","MLST: Sequence type", "Internal UKHSA: Clonal complex (mlst)"]
	krocus = krocus[['Name', "Sample bin", "krocus: Scheme","Sequence type", "Internal UKHSA: Clonal complex"]]
	krocus.columns = ['Name', 'Bin', "krocus: Scheme", "Krocus: Sequence type", "Internal UKHSA: Clonal complex (Krocus)"]

	cc_merged = pd.merge(mlst, krocus,on=['Name', 'Bin'], how = 'outer')
	cc_merged = cc_merged.fillna("NA")
	ccvals = []
	stvals = []
	scheme = []

	for index,row in cc_merged.iterrows():
		if row['Internal UKHSA: Clonal complex (mlst)'] == row['Internal UKHSA: Clonal complex (Krocus)']:
			ccvals.append(row['Internal UKHSA: Clonal complex (mlst)'])
		elif row['Internal UKHSA: Clonal complex (mlst)'] == "NA" and row['Internal UKHSA: Clonal complex (Krocus)'] != "NA":
			ccvals.append(row['Internal UKHSA: Clonal complex (Krocus)'])
		elif row['Internal UKHSA: Clonal complex (mlst)'] != "NA" and row['Internal UKHSA: Clonal complex (Krocus)'] == "NA":
			ccvals.append(row['Internal UKHSA: Clonal complex (mlst)'])
		else:
			ccvals.append("Unknown")

		if row['MLST: Sequence type'] == row['Krocus: Sequence type']:
			stvals.append(row['MLST: Sequence type'])
			scheme.append(row['mlst: Scheme'])
		elif row['MLST: Sequence type'] == "-" and row['Krocus: Sequence type'] == "NA":
			stvals.append("Unknown")
			scheme.append(row['mlst: Scheme'])
		elif row['MLST: Sequence type'] == "-" and row['Krocus: Sequence type'] != "ND":
			stvals.append(row['Krocus: Sequence type'])
			scheme.append(row['krocus: Scheme'])
		elif row['MLST: Sequence type'] != "-" and row['Krocus: Sequence type'] == "ND":
			stvals.append(row['MLST: Sequence type'])
			scheme.append(row['mlst: Scheme'])
		else:
			stvals.append("Unknown")
			scheme.append(row['mlst: Scheme'])

	cc_merged['Clonal complex'] = ccvals
	cc_merged['Sequence type'] = stvals
	cc_merged['Scheme'] = scheme

	cc_merged = cc_merged[['Name', 'Bin', 'Scheme', 'Sequence type', 'Clonal complex']]
	cc_merged.replace("NA","Unknown", inplace=True, regex=False)

	cc_merged['Bin'] = cc_merged['Bin'].str.replace('bin_', '')
	cc_merged.columns = ['Name', 'Bin','scheme', 'st', 'cc']

	cc_merged['st'] = cc_merged['st'].apply(lambda x: str(x).replace('.0',''))
	cc_merged = cc_merged.sort_values('Bin')

	cc_merged_data = cc_merged.to_dict(orient='records')

	return cc_merged_data


def merge_typing(genefinding_filenames, bin_summary_clean):
	typing_merged = pd.concat([pd.read_csv(f, sep='\t') for f in genefinding_filenames ])

	typing_merged = typing_merged[['Bin','Genes identified']]
	typing_merged.columns = ['Bin', 'gi' ]
	typing_merged['Bin'] = typing_merged['Bin'].str.replace('bin_', '')
	typing_merged['gi'] = typing_merged['gi'].str.removesuffix(',')
	typing_merged = typing_merged.sort_values('Bin').fillna("")

	bin_summary_clean = bin_summary_clean[['Bin','itc']]
	typing_merged['Bin'] = typing_merged['Bin'].str.replace('bin_', '')
	typing_merged2 = pd.merge(bin_summary_clean, typing_merged, on=['Bin'], how = 'right')

	typing_summary_data = typing_merged2.to_dict(orient='records')

	return typing_summary_data


def process_amr_typing(hamr_filenames, tool):
	hamr = pd.concat([pd.read_csv(f, sep='\t') for f in hamr_filenames ])

	if (hamr['analysis_software_name'] == tool).any().any():

		if tool == "resfinder":
			hamr1 = hamr[hamr['analysis_software_name'] == tool]
			hamr2 = hamr[hamr['analysis_software_name'] == "pointfinder"]
			hamr = pd.concat([hamr1, hamr2])
		elif tool == "all":
			pass
		else:
			hamr = hamr[hamr['analysis_software_name'] == tool]

		hamr['analysis_software_name'] = hamr['analysis_software_name'].str.replace("resfinder","ResFinder").replace("amrfinderplus","AMRFinderPlus").replace("abricate", "ABRicate").replace("rgi","RGI").replace("pointfinder", "PointFinder")
		hamr['Name'] = hamr['input_file_name'].str.split('.').str[0]
		hamr['Bin'] = hamr['input_file_name'].str.split('.').str[1]
		hamr = hamr[['Name','Bin','gene_symbol','reference_database_name','reference_database_version','analysis_software_name','analysis_software_version','drug_class']]

		merged_amr = hamr.groupby(['Name','Bin', 'reference_database_name', 'reference_database_version', 'analysis_software_name', 'analysis_software_version']).agg(lambda x: ', '.join(x)).reset_index()

		merged_amr['drug_class'] = merged_amr['drug_class'].apply(remove_duplicated_strings)
		merged_amr['gene_symbol'] = merged_amr['gene_symbol'].apply(remove_duplicated_strings)
		merged_amr['Bin'] = merged_amr['Bin'].str.replace('bin_', '')

		amr_summary_data = merged_amr.to_dict(orient='records')

	else:
		amr_summary_data = "None"

	return amr_summary_data


#
def process_plasmidfinder(plasmidfinder_filenames):
	plasmidfinder_summary = pd.concat([pd.read_csv(f, sep='\t') for f in plasmidfinder_filenames ])
	plasmidfinder_summary['Bin'] = plasmidfinder_summary['Contig'].str.split('.').str[1].str.replace('bin_', '')
	plasmidfinder_summary['Note'] = plasmidfinder_summary['Note'].fillna("")
	plasmidfinder_summary = plasmidfinder_summary.sort_values('Bin')

	plasmidfinder_summary = plasmidfinder_summary[['Bin', 'Database', 'Plasmid', 'Identity', 'Note']]
	plasmidfinder_summary_data = plasmidfinder_summary.to_dict(orient='records')

	return plasmidfinder_summary_data


#
def process_ecoli(ecoli_filenames):
	ecoli_summary = pd.concat([pd.read_csv(f, sep='\t') for f in ecoli_filenames ])
	ecoli_summary['Bin'] = ecoli_summary['Bin'].str.replace('bin_', '')

	ecoli_summary = ecoli_summary[['Bin', 'Inferred classification', 'Consensus serotype', 'Consensus ipaH', 'Shigatyper: ipaB', 'Shigeifinder: Virulence_plasmid', 'STECFinder: stx type']]
	ecoli_summary.columns = ['Bin', 'inf_c', 'con_sero', 'con_ipah', 'ipab', 'shig_vp','stx']

	ecoli_summary_data = ecoli_summary.to_dict(orient='records')

	return ecoli_summary_data


#
def process_salmonella(salmonella_filenames):
	salmonella_summary = pd.concat([pd.read_csv(f, sep='\t') for f in salmonella_filenames ])
	salmonella_summary['Bin'] = salmonella_summary['Bin'].str.replace('bin_', '')
	salmonella_summary = salmonella_summary[['Bin', 'Seqsero2: Predicted identification', 'SeqSero2: Predicted serotype', 'SeqSero2: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))', 'SISTR: Predicted antigenic profile (O antigen:H1 antigen (fliC): H2 antigen (fljB))','SISTR: Serogroup', 'SISTR: Serovar antigen']]
	salmonella_summary.columns = ['Bin', 'seqs_pi', 'seqs_ps', 'seqs_pap', 'sistr_pap', 'sistr_sv', 'sistr_sa']
	salmonella_summary['sistr_sa'] = salmonella_summary['sistr_sa'].str.replace('|', ' or ')

	salmonella_summary_data = salmonella_summary.to_dict(orient='records')

	return salmonella_summary_data


#
def process_lmonocytogenes(lmono_filenames):
	lmono_summary = pd.concat([pd.read_csv(f, sep='\t') for f in lmono_filenames ])
	lmono_summary['Bin'] = lmono_summary['ID'].str.split("/").str[-1].str.split(".").str[2].str.replace('bin_', '')
	lmono_summary = lmono_summary[['Bin','SEROTYPE', 'PRS', 'LMO0737', 'LMO1118', 'ORF2110', 'ORF2819', 'COMMENT']]
	lmono_summary = lmono_summary.sort_values('Bin')
	lmono_summary = lmono_summary.replace('FULL', 'Full').replace('NONE', "None").replace("PARTIAL","Partial")

	lmono_summary['COMMENT'] = lmono_summary['COMMENT'].fillna("")
	lmono_summary_data = lmono_summary.to_dict(orient='records')

	return lmono_summary_data

#
def remove_duplicated_strings(s):
	items = s.split(', ')
	unique_items = sorted(set(items), key=items.index)

	return ', '.join(unique_items)


def merge_ccvals(sample_data, readqc_data, bracken_targets_data, contig_data, bin_data, bin_summary_data, run_meta_data, cc_merged, typing_merged, merged_amr, bracken_all_data, ecoli_summary_data, plasmidfinder_summary_data, salmonella_summary_data, lmonocytogenes_summary_data, tax_tool):
	context = {}
	context.update(run_meta_data)
	context.update(sample_data)

	if readqc_data != "None":
		context.update(readqc_data)

	if contig_data != "None":
		context.update(contig_data)

	if bin_data != "None":
		context.update(bin_data)

	if bin_summary_data != "None":
		context['bin_summary'] = bin_summary_data

	if bracken_targets_data != "None":
		context['taxonomic_summary'] = bracken_targets_data
		tax_tool_target = {'tax_tool_target' : tax_tool}
		context.update(tax_tool_target)

	if bracken_all_data != "None":
		context['taxonomic_summary_all'] = bracken_all_data
		tax_tool_all = {'tax_tool_all' : tax_tool}
		context.update(tax_tool_all)

	if cc_merged != "None":
		context['mlst_summary'] = cc_merged

	if typing_merged != "None":
		context['typing_summary'] = typing_merged

	if merged_amr != "None":
		context['amr_summary'] = merged_amr

	if ecoli_summary_data != "None":
		context['ecoli_summary'] = ecoli_summary_data

	if plasmidfinder_summary_data != "None":
		context['plasmid_summary'] = plasmidfinder_summary_data

	if salmonella_summary_data != "None":
		context['salmonella_summary'] = salmonella_summary_data

	if lmonocytogenes_summary_data != "None":
		context['lmonocytogenes_summary'] = lmonocytogenes_summary_data

	return context


def render_template(context, template, output):
	template_dir = os.getcwd()
	env = Environment(loader=FileSystemLoader(template_dir))
	template = env.get_template(template)

	html_output = template.render(context)

	html_output = html_output.replace("<em>Unclassified</em>","Unclassified")
	html_output = html_output.replace("<em>Unclassified Bacteria</em>","Unclassified Bacteria")

	with open(output + '.summary_report.html', 'w') as f:
		f.write(html_output)


# Parse arguments from the command line.
def parse_args():
	description = 'Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--sample_id', required=True, help="Sample ID")
	parser.add_argument('--run_id', required=True, help="Run ID")
	parser.add_argument('--barcode', required=True, help="Barcode")
	parser.add_argument('--sample_type', required=True, help="Sample type")
	parser.add_argument('--readqc', required=False, help="Read QC stats")
	parser.add_argument('--hostreads', required=False, help="List of reads flagged as being derived from host")
	parser.add_argument('--kraken2_bracken_targets', required=False, help="Bracken output, subset for target species")
	parser.add_argument('--kraken2_bracken_all', required=False, help="Bracken output")
	parser.add_argument('--centrifuger_bracken_targets', required=False, help="Bracken output, subset for target species")
	parser.add_argument('--centrifuger_bracken_all', required=False, help="Bracken output")

	parser.add_argument('--tax_mode', required=False, help="Taxonomic assignment results to report")
	parser.add_argument('--min_read_prop', required=False, help="Min proportion of reads required to retain taxonomic hit", default=0)
	parser.add_argument('--min_read_count', required=False, help="Min count of reads required to retain taxonomic hit", default=0)
	parser.add_argument('--contig_summary', required=False, help="Contig summary file")
	parser.add_argument('--bin_summary', required=False, help='Bin summary file')
	parser.add_argument('--krocus', required=False, nargs = '*', help='Krocus_ccvals')
	parser.add_argument('--mlst', required=False, nargs = '*', help='MLST ccvals')
	parser.add_argument('--genefinding', required=False, nargs = '*', help='Genefinding results')
	parser.add_argument('--plasmidfinder', required=False, nargs = '*', help='Plasmidfinder results')
	parser.add_argument('--ecoli_typing', required=False, nargs = '*', help='E. coli typing results')
	parser.add_argument('--salmonella_typing', required=False, nargs = '*', help='Salmonella typing results')
	parser.add_argument('--lmonocytogenes_typing', required=False, nargs = '*', help='L. monocytogenes typing results')
	parser.add_argument('--hamronization_summary', required=False, nargs = '*', help = 'hAMRonization combined report')
	parser.add_argument('--amr_tool', required=False, help='Tool to report AMR detection')
	parser.add_argument('--pipeline_version', required=True, help='Pipeline version')
	parser.add_argument('--logo', required=False, help="Logo")
	parser.add_argument('--report_template', required=True, help="HTML template")
	parser.add_argument('--output', required=True, help="Output HTML filename prefix")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()
	run_meta_data = process_pipeline_metadata(args.pipeline_version)
	sample_data = process_metadata(args.sample_id, args.run_id, args.barcode, args.sample_type, args.logo)

	if args.readqc and args.hostreads:
		readqc_data = process_readqc(args.readqc, args.hostreads)
	else:
		readqc_data = "None"

	if args.kraken2_bracken_targets and args.min_read_count and args.min_read_prop and args.tax_mode == "Kraken2":
		bracken_targets_data = process_bracken_targets(args.kraken2_bracken_targets, args.min_read_count, args.min_read_prop)
	else:
		bracken_targets_data = "None"

	if args.kraken2_bracken_all and args.min_read_count and args.min_read_prop and args.tax_mode == "Kraken2":
		bracken_all_data = process_bracken_all(args.kraken2_bracken_all, args.min_read_count, args.min_read_prop)
	else:
		bracken_all_data = "None"

	if args.centrifuger_bracken_targets and args.min_read_count and args.min_read_prop and args.tax_mode == "Centrifuger":
		bracken_targets_data = process_bracken_targets(args.centrifuger_bracken_targets, args.min_read_count, args.min_read_prop)
	if args.centrifuger_bracken_all and args.min_read_count and args.min_read_prop and args.tax_mode == "Centrifuger":
		bracken_all_data = process_bracken_all(args.centrifuger_bracken_all, args.min_read_count, args.min_read_prop)

	if args.bin_summary:
		bin_data, bin_summary = process_bin_summary(args.bin_summary)
		bin_summary_data, bin_summary_clean = process_bin_typing(bin_summary)
	else:
		bin_data = "None"
		bin_summary_data = "None"
		bin_summary_clean = "None"

	if args.contig_summary:
		contig_data = process_contig_summary(args.contig_summary)
	else:
		contig_data = "None"

	if args.krocus and args.mlst:
		cc_merged = process_mlst(args.krocus, args.mlst)
	elif args.krocus and not args.mlst:
		cc_merged = process_krocus_indv(args.krocus)
	elif args.mlst and not args.krocus:
		cc_merged = process_mlst_indv(args.mlst)
	else:
		cc_merged = "None"

	if args.genefinding:
		typing_summary_data = merge_typing(args.genefinding, bin_summary_clean)
	else:
		typing_summary_data = "None"

	if args.hamronization_summary and args.amr_tool:
		amr_summary_data = process_amr_typing(args.hamronization_summary, args.amr_tool)
	else:
		amr_summary_data = "None"

	if args.plasmidfinder:
		plasmidfinder_summary_data = process_plasmidfinder(args.plasmidfinder)
	else:
		plasmidfinder_summary_data = "None"

	if args.ecoli_typing:
		ecoli_summary_data = process_ecoli(args.ecoli_typing)
	else:
		ecoli_summary_data = "None"

	if args.salmonella_typing:
		salmonella_summary_data = process_salmonella(args.salmonella_typing)
	else:
		salmonella_summary_data = "None"

	if args.lmonocytogenes_typing:
		lmonocytogenes_summary_data = process_lmonocytogenes(args.lmonocytogenes_typing)
	else:
		lmonocytogenes_summary_data = "None"

	context = merge_ccvals(sample_data, readqc_data, bracken_targets_data, contig_data, bin_data, bin_summary_data, run_meta_data, cc_merged, typing_summary_data, amr_summary_data, bracken_all_data, ecoli_summary_data, plasmidfinder_summary_data, salmonella_summary_data, lmonocytogenes_summary_data, args.tax_mode)
	render_template(context, args.report_template, args.output)


if __name__ == "__main__":
	main()
