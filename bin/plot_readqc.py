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
import math
import seaborn as sns
import matplotlib.pyplot as plt
from pretty_html_table import build_table
import base64
import matplotlib.image as mpimg
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
			'logo' : base64.b64encode(image_file.read()).decode('utf-8')
		}

	return sample_data


# Process and reformat Nanoplot output table (nanostats.txt)
def process_nanostats(nanostats_in, stage):
	df_ns = pd.read_csv(nanostats_in, sep='|', header=None)
	df_ns = df_ns[df_ns[0].str.contains(":")]
	df_ns = df_ns[~df_ns[0].str.contains("Top")]
	df_ns = df_ns[~df_ns[0].str.contains("percentage")]
	df_ns = df_ns[~df_ns[0].str.contains("General")]
	df_ns_1 = pd.concat([df_ns[0].str.split(':', expand=True).add_prefix('SubColumn')], axis=1)

	df_ns_2 = pd.concat([df_ns_1['SubColumn0'], df_ns_1['SubColumn1'].str.split('\t', expand=True).add_prefix('SubColumn')], axis=1)
	df_ns_3 = df_ns_2.head(8)
	replacement = {
		"Mean read length": "Mean read length (bp)",
		"Median read length": "Median read length (bp)",
		"Read length N50": u"Read N\u2085\u2080 (bp)",
		"STDEV read length": "Standard deviation of read length (bp)",
		"Total bases": "Total bases (Mbp)"}

	df_ns_4 = df_ns_3['SubColumn0'].replace(replacement)
	df_ns_4.columns = ['', stage]
	df_ns_4.style.hide(axis='index')

	return(df_ns_4)


# Merge pre- and post-QC nanostats results
def merge_nanostats(ns_preqc, ns_postqc, output):
	ns_merged = pd.merge(ns_preqc, ns_postqc, on="", how="inner").replace(",","", regex=True)
	ns_merged = ns_merged.loc[[4,0,2,6,5,1,3,7]]
	ns_merged['Change (%)'] = (((pd.to_numeric(ns_merged['Post-QC'])-pd.to_numeric(ns_merged['Pre-QC']))/pd.to_numeric(ns_merged['Pre-QC']))*100).round(2)

	ns_merged2 = ns_merged
	ns_merged2['Pre-QC'] = ns_merged2['Pre-QC'].str.strip()
	ns_merged2['Post-QC'] = ns_merged2['Post-QC'].str.strip()

	ns_merged2['Pre-QC'] = pd.to_numeric(ns_merged2['Pre-QC'])
	ns_merged2['Post-QC'] = pd.to_numeric(ns_merged2['Post-QC'])

	ns_merged2.loc[[7], "Pre-QC"] = round((ns_merged2.loc[[7], "Pre-QC"]/1000000),2)
	ns_merged2.loc[[7], "Post-QC"] = round((ns_merged2.loc[[7], "Post-QC"]/1000000),2)

	ns_merged2.loc[[1], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[1])
	ns_merged2.loc[[1], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[1])
	ns_merged2.loc[[2], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[2])
	ns_merged2.loc[[2], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[2])
	ns_merged2.loc[[3], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[3])
	ns_merged2.loc[[3], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[3])
	ns_merged2.loc[[4], "Pre-QC"] = '{:.0f}'.format(ns_merged2['Pre-QC'].loc[4])
	ns_merged2.loc[[4], "Post-QC"] = '{:.0f}'.format(ns_merged2['Post-QC'].loc[4])
	ns_merged2.loc[[5], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[5])
	ns_merged2.loc[[5], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[5])
	ns_merged2.loc[[6], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[6])
	ns_merged2.loc[[6], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[6])
	ns_merged2.loc[[0], "Pre-QC"] = '{:.2f}'.format(ns_merged2['Pre-QC'].loc[0])
	ns_merged2.loc[[0], "Post-QC"] = '{:.2f}'.format(ns_merged2['Post-QC'].loc[0])

	ns_merged2.to_csv(output + ".nanostats.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')

	return(ns_merged2)


# Process Seqtk fqchk results for the starts and ends of reads
def process_seqtk_fqchk(fqchk_in_fw, fqchk_in_rv):
	fqchk_fw = pd.read_csv(fqchk_in_fw, sep='\t', header=None)
	fqchk_fw['Orientation'] = "FW"

	fqchk_rv = pd.read_csv(fqchk_in_rv, sep='\t', header=None)
	fqchk_rv['Orientation'] = "RV"

	fqchk_merged = pd.concat([fqchk_rv, fqchk_fw])
	fqchk_merged.columns = ["POS","#bases","%A","%C","%G","%T","%N","avgQ","errQ","%low","%high","Orientation"]
	fqchk_merged = fqchk_merged[["POS", "%A","%C","%T","%G","avgQ","Orientation"]]

	return(fqchk_merged)


# Plot Seqtk fqchk results for both pre- and post-QC
def plot_seqtk_fqchk(fqchk_merged, stage):
	maxq = math.ceil((max(fqchk_merged["avgQ"])*1.2)/5)*5
	fqchk_fw = fqchk_merged[fqchk_merged['Orientation'] == "FW"]
	fqchk_rv = fqchk_merged[fqchk_merged['Orientation'] == "RV"]

	fw_fig, fw_ax1 = plt.subplots()
	fw_ax1.bar(fqchk_fw['POS'], fqchk_fw['%G'], width=1, bottom=fqchk_fw['%A']+fqchk_fw['%T']+fqchk_fw['%C'],color='#AE93BE')
	fw_ax1.bar(fqchk_fw['POS'], fqchk_fw['%C'], width=1, bottom=fqchk_fw['%A']+fqchk_fw['%T'], color='#B4DAE5')
	fw_ax1.bar(fqchk_fw['POS'], fqchk_fw['%T'], width=1, bottom=fqchk_fw['%A'], color='#F0D77B')
	fw_ax1.bar(fqchk_fw['POS'], fqchk_fw['%A'], width=1, color='#bab8b8')
	fw_ax1.set_ylim([0, 100])
	fw_ax1.set_xlim([0, 300])
	fw_ax1.xaxis.set_ticks(np.arange(0, 301, 100))
	fw_ax1.set_ylabel('Nucleotide composition (%)', weight='bold')
	fw_ax1.set_xlabel("5' position (bp)", weight='bold')
	fw_ax1.legend(['G (%)','C (%)','T (%)','A (%)'],loc='lower right')

	fw_ax2 = fw_ax1.twinx()
	fw_ax2.plot(fqchk_fw['POS'], fqchk_fw["avgQ"], color="red", linewidth=1.6)
	fw_ax2.set_ylabel('Mean PHRED score', color="red", weight='bold')
	fw_ax2.tick_params(axis='y', colors="red")
	fw_ax2.spines['right'].set_color("red")
	fw_ax2.yaxis.set_ticks(np.linspace(0,maxq,4))
	fw_fig.tight_layout()

	if stage == "Pre-QC":
		fw_fig.text(-0.05, 1, "A",size=15, weight='bold')
	if stage == "Post-QC":
		fw_fig.text(-0.05, 1, "C",size=15, weight='bold')
	fw_fig.savefig('fw_'+stage+'.png', bbox_inches="tight")

	rv_fig, rv_ax1 = plt.subplots()
	rv_ax1.bar(fqchk_rv['POS'], fqchk_rv['%G'], width=1, bottom=fqchk_rv['%A']+fqchk_rv['%T']+fqchk_rv['%C'],color='#AE93BE')
	rv_ax1.bar(fqchk_rv['POS'], fqchk_rv['%C'], width=1, bottom=fqchk_rv['%A']+fqchk_rv['%T'], color='#B4DAE5')
	rv_ax1.bar(fqchk_rv['POS'], fqchk_rv['%T'], width=1, bottom=fqchk_rv['%A'], color='#F0D77B')
	rv_ax1.bar(fqchk_rv['POS'], fqchk_rv['%A'], width=1, color='#bab8b8')
	rv_ax1.set_ylim([0, 100])
	rv_ax1.set_xlim([0, 300])
	rv_ax1.xaxis.set_ticks(np.arange(0, 301, 100))
	rv_ax1.set_ylabel('Nucleotide composition (%)', weight='bold')
	rv_ax1.set_xlabel("3' position (bp)", weight='bold')

	rv_ax2 = rv_ax1.twinx()
	rv_ax2.plot(fqchk_rv['POS'], fqchk_rv["avgQ"], color="red", linewidth=1.6)
	rv_ax2.set_ylabel('Mean PHRED score', color="red", weight='bold')
	rv_ax2.tick_params(axis='y', colors="red")
	rv_ax2.spines['right'].set_color("red")
	rv_ax2.yaxis.set_ticks(np.linspace(0,maxq,4))
	rv_ax1.invert_xaxis()
	rv_fig.tight_layout()

	if stage == "Pre-QC":
		rv_fig.text(-0.05, 1, "B",size=15, weight='bold')
	if stage == "Post-QC":
		rv_fig.text(-0.05, 1, "D",size=15, weight='bold')

	rv_fig.savefig('rv_'+stage+'.png', bbox_inches="tight")


	f, axarr = plt.subplots(1,2, figsize=(14, 14))
	axarr[0].imshow(mpimg.imread('fw_'+stage+'.png'))
	axarr[1].imshow(mpimg.imread('rv_'+stage+'.png'))
	[ax.set_axis_off() for ax in axarr.ravel()]
	plt.tight_layout()

	buffer = BytesIO()
	plt.savefig(buffer, bbox_inches="tight", format="png")
	buffer.seek(0)
	image_base64 = base64.b64encode(buffer.read()).decode('utf-8')
	buffer.close()

	stage2 = stage.replace("-QC","")

	fig_data = {
		'image_base64_'+stage2 : image_base64,
	}

	return(fig_data)


# Combine read QC statistics
def process_readqc(read_qc_stats):
	read_qc_stats.columns = ["key", "pre", "post", "change"]
	readqc_data = read_qc_stats.to_dict(orient='records')

	return readqc_data


# Process raw Nanoplot data (NanoPlot-data.tsv.gz)
def process_nanodata(fn_pre, fn_post):
	nd_pre = pd.read_csv(fn_pre, compression='gzip', sep='\t')
	nd_pre['Stage'] = "Pre-QC"

	nd_post = pd.read_csv(fn_post, compression='gzip', sep='\t')
	nd_post['Stage'] = "Post-QC"

	nd_merged = pd.concat([nd_pre, nd_post], ignore_index=True, sort=False)
	nd_merged["klen"] = nd_merged["lengths"]/1000

	return(nd_merged)


# Plot raw Nanoplot data and merge subplots
def plot_nanodata(nd_merged):
	re_nd_pre = nd_merged.loc[nd_merged['Stage']== "Pre-QC"]
	re_nd_post = nd_merged.loc[nd_merged['Stage']== "Post-QC"]
	roundedqual = math.ceil(nd_merged['quals'].max()/5)*5
	roundedlen = (math.ceil(nd_merged['klen'].max()/50)*50)
	gfg_pre = sns.jointplot(
		data=re_nd_pre,
		x="klen",
		alpha=.25,
		color="r",
		s=5,
		y="quals",
		kind='scatter',
		marginal_kws=dict(bins=100),
		xlim = (0,roundedlen), ylim = (0,roundedqual))
	gfg_pre.set_axis_labels("Read length (kbp)","Mean read quality", fontweight='bold')
	gfg_pre.figure.tight_layout()

	gfg_post = sns.jointplot(
		data=re_nd_post,
		x="klen",
		alpha=.25,
		s=5,
		y="quals",
		kind='scatter',
		marginal_kws=dict(bins=100),
		xlim = (0,roundedlen), ylim = (0,roundedqual))
	gfg_post.set_axis_labels("Read length (kbp)","Mean read quality", fontweight='bold')
	gfg_post.figure.tight_layout()

	gfg_pre.savefig("gfg_pre.png")
	gfg_post.savefig("gfg_post.png")

	f, axarr = plt.subplots(1,2, figsize=(14, 14))
	axarr[0].imshow(mpimg.imread('gfg_pre.png'))
	axarr[1].imshow(mpimg.imread('gfg_post.png'))
	[ax.set_axis_off() for ax in axarr.ravel()]
	plt.tight_layout()

	buffer = BytesIO()
	plt.savefig(buffer, bbox_inches="tight", format="png")
	buffer.seek(0)
	image_base64 = base64.b64encode(buffer.read()).decode('utf-8')
	buffer.close()

	ns_rawfig_data = {
		'ns_raw_image_base64' : image_base64,
	}

	return(ns_rawfig_data)


def merge_ccvals(sample_data, readqc_data, pre_fig_data, post_fig_data, ns_rawfig_data):
	context = {}
	context.update(sample_data)

	if pre_fig_data != "None":
		context.update(pre_fig_data)
	if post_fig_data != "None":
		context.update(post_fig_data)
	if ns_rawfig_data != "None":
		context.update(ns_rawfig_data)
	if readqc_data != "None":
		context['readqc_summary'] = readqc_data

	return context


def render_template(context, template, output):
	template_dir = os.getcwd()
	env = Environment(loader=FileSystemLoader(template_dir))
	template = env.get_template(template)

	html_output = template.render(context)

	with open(output + '.readqc_report.html', 'w') as f:
		f.write(html_output)


# Parse arguments from the command line.
def parse_args():
	description = 'Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description="description")
	parser.add_argument('--sample_id', required=True, help="Sample ID")
	parser.add_argument('--run_id', required=True, help="Run ID")
	parser.add_argument('--barcode', required=True, help="Barcode")
	parser.add_argument('--sample_type', required=True, help="Sample type")
	parser.add_argument('--nanostats_pre', type=str, help='Nanostats file (pre-QC)')
	parser.add_argument('--nanostats_post', type=str, help='Nanostats file (post-QC)')
	parser.add_argument('--nucl_comp_pre_fw', type=str, help="Nucleotide compositions, 5' (pre-QC)")
	parser.add_argument('--nucl_comp_pre_rv', type=str, help="Nucleotide compositions, 3' (pre-QC)")
	parser.add_argument('--nucl_comp_post_fw', type=str, help="Nucleotide compositions, 5' (post-QC)")
	parser.add_argument('--nucl_comp_post_rv', type=str, help="Nucleotide compositions, 3' (post-QC)")
	parser.add_argument('--nanoplot_raw_pre', type=str, help='Nanoplot raw results (pre-QC)')
	parser.add_argument('--nanoplot_raw_post', type=str, help='Nanoplot raw results (post-QC)')
	parser.add_argument('--logo', required=False, help="Logo")
	parser.add_argument('--report_template', required=True, help="HTML template")
	parser.add_argument('--output', required=True, help="Output HTML filename prefix")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()
	sample_data = process_metadata(args.sample_id, args.run_id, args.barcode, args.sample_type, args.logo)

	if args.nanostats_pre and args.nanostats_post:
		ns_preqc = process_nanostats(args.nanostats_pre, "Pre-QC")
		ns_postqc = process_nanostats(args.nanostats_post, "Post-QC")
		ns_merged = merge_nanostats(ns_preqc, ns_postqc, args.output)
		readqc_data = process_readqc(ns_merged)
	else:
		readqc_data = "None"

	if args.nucl_comp_pre_fw and args.nucl_comp_pre_rv:
		fqchk_merged_pre = process_seqtk_fqchk(args.nucl_comp_pre_fw, args.nucl_comp_pre_rv)
		pre_fig_data = plot_seqtk_fqchk(fqchk_merged_pre, "Pre-QC")
		fqchk_merged_post = process_seqtk_fqchk(args.nucl_comp_post_fw, args.nucl_comp_post_rv)
		post_fig_data = plot_seqtk_fqchk(fqchk_merged_post, "Post-QC")
	else:
		pre_fig_data = "None"
		post_fig_data = "None"

	if args.nanoplot_raw_pre and args.nanoplot_raw_post:
		df4 = process_nanodata(args.nanoplot_raw_pre, args.nanoplot_raw_post)
		ns_rawfig_data = plot_nanodata(df4)
	else:
		ns_rawfig_data = "None"

	context = merge_ccvals(sample_data, readqc_data, pre_fig_data, post_fig_data, ns_rawfig_data)
	render_template(context, args.report_template, args.output)


if __name__ == "__main__":
	main()

