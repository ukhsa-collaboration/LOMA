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

from io import BytesIO


###### Functions


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


# Summarize binning statistics
def process_bin_summary(bin_summary_ccvals):
	bin_summary = pd.read_csv(bin_summary_ccvals, sep='\t')
	bin_summary['Bin'] = bin_summary['Bin'].str.replace('bin_', '')
	bin_summary = bin_summary[['Bin', 'Inferred taxononomic classification','Inferred classification rank']]
	bin_summary['Bin'] = bin_summary['Bin'].str.replace('bin_', '')
	bin_summary.columns = ['Bin', 'itc', 'icr']

	return(bin_summary)



#
def process_amr_typing(hamr_filenames, tool, bin_summary):
	hamr = pd.concat([pd.read_csv(f, sep='\t') for f in hamr_filenames ])

	if (hamr['analysis_software_name'] == tool).any().any():
		if tool == "resfinder":
			hamr1 = hamr[hamr['analysis_software_name'] == tool]
			hamr2 = hamr[hamr['analysis_software_name'] == "pointfinder"]
			hamr = pd.concat([hamr1, hamr2])
		else:
			hamr = hamr[hamr['analysis_software_name'] == tool]

		hamr['Name'] = hamr['input_file_name'].str.split('.').str[0]
		hamr['Bin'] = hamr['input_file_name'].str.split('.').str[1]
		hamr['Bin'] = hamr['Bin'].str.replace('bin_', '')
		hamr['reference_database_name'] = hamr['reference_database_name'].str.replace("ncbi","NCBI")
		hamr['analysis_software_name'] = hamr['analysis_software_name'].str.replace("resfinder","ResFinder").replace("amrfinderplus","AMRFinderPlus").replace("abricate", "ABRicate").replace("rgi","RGI").replace("pointfinder","PointFinder")
		hamr['genetic_variation_type'] = hamr['genetic_variation_type'].str.replace('protein_variant_detected', 'Protein variant').str.replace('gene_presence_detected', 'Gene presence')
		hamr['drug_class'] = hamr['drug_class'].str.lower().str.title()

		hamr2 = hamr[['Name','Bin','gene_symbol','gene_name','reference_database_name','reference_database_version','analysis_software_name','analysis_software_version','genetic_variation_type','drug_class','amino_acid_mutation','resistance_mechanism']]
		hamr2 = hamr2.fillna("")
		hamr = hamr[['Name','Bin','gene_symbol','reference_database_name','reference_database_version','analysis_software_name','analysis_software_version','drug_class']]

		merged_amr = hamr.groupby(['Name','Bin', 'reference_database_name', 'reference_database_version', 'analysis_software_name', 'analysis_software_version']).agg(lambda x: ', '.join(x)).reset_index()

		merged_amr['drug_class'] = merged_amr['drug_class'].apply(remove_duplicated_strings)
		merged_amr['gene_symbol'] = merged_amr['gene_symbol'].apply(remove_duplicated_strings)

		merged_amr = merged_amr.sort_values('Bin')
		hamr2 = hamr2.sort_values('Bin')

		merged_amr2 = pd.merge(bin_summary, merged_amr, on=['Bin'], how = 'right')
		hamr3 = pd.merge(bin_summary, hamr2, on=['Bin'], how = 'right')

		amr_summary = merged_amr2.to_dict(orient='records')
		detailed_summary = hamr3.to_dict(orient='records')

	else:
		amr_summary = "None"
		detailed_summary = "None"

	return(amr_summary, detailed_summary)


#
def remove_duplicated_strings(s):
	items = s.split(', ')
	unique_items = sorted(set(items), key=items.index)

	return ', '.join(unique_items)


def merge_ccvals(sample_data, resfinder_merged_amr, resfinder_detailed_summary, amrfinderplus_merged_amr, amrfinderplus_detailed_summary, abricate_merged_amr, abricate_detailed_summary, rgi_merged_amr, rgi_detailed_summary):
	context = {}
	context.update(sample_data)

	if resfinder_merged_amr != "None":
		context['resfinder_amr_summary'] = resfinder_merged_amr
		context['resfinder_detailed_summary'] = resfinder_detailed_summary

	if amrfinderplus_merged_amr != "None":
		context['amrfinderplus_amr_summary'] = amrfinderplus_merged_amr
		context['amrfinderplus_detailed_summary'] = amrfinderplus_detailed_summary

	if abricate_merged_amr != "None":
		context['abricate_amr_summary'] = abricate_merged_amr
		context['abricate_detailed_summary'] = abricate_detailed_summary

	if rgi_merged_amr != "None":
		context['rgi_amr_summary'] = rgi_merged_amr
		context['rgi_detailed_summary'] = rgi_detailed_summary

	return context


def render_template(context, template, output):
	template_dir = os.getcwd()
	env = Environment(loader=FileSystemLoader(template_dir))
	template = env.get_template(template)

	html_output = template.render(context)

	html_output = html_output.replace("<em>Unclassified</em>","Unclassified")
	html_output = html_output.replace("<em>Unclassified Bacteria</em>","Unclassified Bacteria")

	with open(output + '.amr_report.html', 'w') as f:
		f.write(html_output)


# Parse arguments from the command line.
def parse_args():
	description = 'Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--sample_id', required=True, help="Sample ID")
	parser.add_argument('--run_id', required=True, help="Run ID")
	parser.add_argument('--barcode', required=True, help="Barcode")
	parser.add_argument('--sample_type', required=True, help="Sample type")
	parser.add_argument('--hamronization_summary', required=False, nargs = '*', help = 'hAMRonization combined report')
	parser.add_argument('--amr_tool', required=False, help='Tool to report AMR detection')
	parser.add_argument('--bin_summary', required=False, help='Bin summary file')
	parser.add_argument('--logo', required=False, help="Logo")
	parser.add_argument('--report_template', required=True, help="HTML template")
	parser.add_argument('--output', required=True, help="Output HTML filename prefix")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()
	sample_data = process_metadata(args.sample_id, args.run_id, args.barcode, args.sample_type, args.logo)

	if args.bin_summary:
		bin_summary = process_bin_summary(args.bin_summary)

	if args.hamronization_summary:
		resfinder_amr_summary, resfinder_detailed_summary = process_amr_typing(args.hamronization_summary, "resfinder", bin_summary)
		amrfinderplus_amr_summary, amrfinderplus_detailed_summary = process_amr_typing(args.hamronization_summary, "amrfinderplus", bin_summary)
		abricate_amr_summary, abricate_detailed_summary = process_amr_typing(args.hamronization_summary, "abricate", bin_summary)
		rgi_amr_summary, rgi_detailed_summary = process_amr_typing(args.hamronization_summary, "rgi", bin_summary)
	else:
		resfinder_amr_summary = "None"
		resfinder_detailed_summary = "None"
		amrfinderplus_amr_summary = "None"
		amrfinderplus_detailed_summary = "None"
		abricate_amr_summary = "None"
		abricate_detailed_summary = "None"
		rgi_amr_summary = "None"
		rgi_detailed_summary = "None"

	context = merge_ccvals(sample_data, resfinder_amr_summary, resfinder_detailed_summary, amrfinderplus_amr_summary, amrfinderplus_detailed_summary, abricate_amr_summary, abricate_detailed_summary, rgi_amr_summary, rgi_detailed_summary)
	render_template(context, args.report_template, args.output)


if __name__ == "__main__":
	main()
