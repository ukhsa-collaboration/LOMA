#!/usr/bin/env python

# Plots individual and (optionally) merged input from read-based taxonomic classifiers.
# Designed to run with Kraken2 (Bracken and/or Taxpasta), Centrifuger (Bracken and/or Taxpasta) and Sylph results.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports

import pandas as pd
import numpy as np
import plotly.express as px
import os
import argparse
from jinja2 import Environment, FileSystemLoader
import base64
from io import BytesIO
import plotly.io as pio


###### Functions

plot_B_colors = ["#ff4a24","#19D3F3","#F58518","#e60019","#006fdf","#FECB52","#636EFA","#80ba5a","#b596ed","#009b3b","#bc80bd","#882255","#117733","#88ccee","#ff0b6e","#00d395","#f2b701","#999933","#E45756","#332288","#cf1c90","#FFA15A","#ccebc5","#80b1d3","#FF9DA6","#002b90","#008695","#11a579","#3969ac","#44aa99","#54A24B","#72B7B2","#7f3c8d","#9D755D","#aa4499","#AB63FA","#B279A2","#b3de69","#B6E880","#bebada","#c44a86","#cc6677","#ddcc77","#de57ff","#e59800","#e68310","#e73f74","#EECA3B","#EF553B","#fb8072","#fccde5","#FF6692","#FF97FF"]
plot_A_colors = ['#ddcc77','#882255','#117733','#332288','#aa4499','#44aa99','#88ccee','#cc6677']


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


#
def process_taxpasta_A(infile_taxp):
	df_taxp = pd.read_csv(infile_taxp, delimiter='\t')

	df_taxp['rank'] = df_taxp['rank'].replace(['kingdom'], 'Kingdom')
	df_taxp['rank'] = df_taxp['rank'].replace(['superkingdom','clade','domain'], 'Domain')
	df_taxp['rank'] = df_taxp['rank'].replace(['tribe','family'], 'Family')
	df_taxp['rank'] = df_taxp['rank'].replace(['superorder','class'], 'Class')
	df_taxp['rank'] = df_taxp['rank'].replace(['superfamily','parvorder','infraorder','suborder','order'], 'Order')
	df_taxp['rank'] = df_taxp['rank'].replace(['superclass','subphylum','phylum'], 'Phylum')
	df_taxp['rank'] = df_taxp['rank'].replace(['superorder','subclass','class'], 'Class')
	df_taxp['rank'] = df_taxp['rank'].replace(['subfamily','family'], 'Family')
	df_taxp['rank'] = df_taxp['rank'].replace(['species'], 'Species')
	df_taxp['rank'] = df_taxp['rank'].replace(['strain','serotype','serogroup','isolate','forma specialis','subspecies'], 'Subspecies')
	df_taxp['rank'] = df_taxp['rank'].replace(['species subgroup','species group','genus','subgenus'], 'Genus')
	df_taxp['rank'] = df_taxp['rank'].replace(['no rank'], 'Missing rank')
	df_taxp['rank'] = df_taxp['rank'].fillna('Unclassified')

	df_taxp_grouped = df_taxp.groupby('rank').sum('count')
	sorter = ['Unclassified', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species','Subspecies','Missing rank']
	df_taxp_grouped = df_taxp_grouped.reindex(index = sorter)
	taxp_total_reads = df_taxp_grouped['count'].sum()
	df_taxp_grouped['perc'] = ((df_taxp_grouped['count']/taxp_total_reads)*100).round(2)
	df_taxp_grouped['x_value'] = 1
	df_taxp_grouped['rank'] = df_taxp_grouped.index

	return(df_taxp, df_taxp_grouped)


#
def plot_taxpasta_A(df_taxp_grouped, tool, output):
	taxp_fig_A = px.bar(df_taxp_grouped,
				x='count',
				y='x_value',
				color='rank',
				color_discrete_sequence=plot_A_colors,
				orientation='h',
				custom_data=['perc', 'count', 'rank'],
				labels={'count': 'Reads' , 'name': 'Species', 'perc': 'Percentage'})

	taxp_fig_A.update_traces(
		hovertemplate="<br>".join([
			"<b>Taxonomic rank</b>: %{customdata[2]}",
			"<b>Reads</b>: %{customdata[1]} (%{customdata[0]}%)",
		] )+"<extra></extra>")

	taxp_fig_A.for_each_trace(lambda t: t.update(name = '<b>' + t.name +'</b>'))
	for t in taxp_fig_A.data:
		taxp_fig_A.update_traces(marker={"color": np.where(t.y>0, "darkgrey", t.marker.color)},
							name = '<b>Unclassified</b>', 
							selector={"name":"<b>Unclassified</b>"})
		taxp_fig_A.update_traces(marker={"color": np.where(t.y>0, "lightgrey", t.marker.color)},
							name = '<b>Missing rank</b>', 
							selector={"name":"<b>Missing rank</b>"})

	taxp_fig_A.update_xaxes(title_text = '<b> Reads (count) </b>')
	taxp_fig_A.update_yaxes(title_text=None, showticklabels=False)
	taxp_fig_A.update_layout(legend_title_text='<b>Taxonomic rank</b>',
						yaxis=dict(range=[0.6,1]),
						plot_bgcolor='white',
						title_text='<b>'+output+'</b>', title_font=dict(size=20))

	taxp_fig_A_string = pio.to_html(taxp_fig_A, full_html=False)
	taxp_fig_A_data = {
		'taxp_fig_A_'+tool : taxp_fig_A_string,
	}

	return(taxp_fig_A_data)


#
def process_taxpasta_B(df_taxp_grouped):
	filtered_df_taxp = df_taxp_grouped[df_taxp_grouped['rank'] == "Species"]
	taxp_rcount = filtered_df_taxp['count'].sum()
	taxp_top_20 = filtered_df_taxp.nlargest(50,'count')
	taxp_other = filtered_df_taxp[filtered_df_taxp['name'] != 'unclassified']
	taxp_other = filtered_df_taxp[~filtered_df_taxp.index.isin(taxp_top_20.index)]
	taxp_other_total_reads = taxp_other['count'].sum()
	taxp_other = pd.DataFrame({'name': ['Other (grouped)'], 'count': [taxp_other_total_reads]})
	taxp_unclassified = df_taxp_grouped[df_taxp_grouped['name'] == 'unclassified']
	df_taxp_stacked = pd.concat([taxp_top_20, taxp_other, taxp_unclassified])
	df_taxp_stacked['perc'] = ((df_taxp_stacked['count']/taxp_rcount)*100).round(2)
	df_taxp_stacked['x_value'] = 1

	return(df_taxp_stacked)


#
def plot_taxpasta_B(df_taxp_stacked, tool, output):
	taxp_fig_B = px.bar(df_taxp_stacked,
						x='count',
						y='x_value',
						color='name',
						color_discrete_sequence=plot_B_colors,
						orientation='h',
						custom_data=['perc', 'count', 'name'],
						labels={'count': 'Reads' , 'name': 'Species', 'perc': 'Percentage'})

	taxp_fig_B.update_layout(legend_title_text='<b>Species (50 most abundant)</b>')
	taxp_fig_B.for_each_trace(lambda t: t.update(name = '<i><b>' + t.name +'</i></b>'))

	for t in taxp_fig_B.data:
		taxp_fig_B.update_traces(marker={"color": np.where(t.y>0, "lightgrey", t.marker.color)},
								name = '<b>Other (grouped)</b>', 
								selector={"name":"<i><b>Other (grouped)</i></b>"})

	taxp_fig_B.update_traces(
		hovertemplate="<br>".join([
			"<b>Species</b>: <i>%{customdata[2]}</i>",
			"<b>Reads</b>: %{customdata[1]} (%{customdata[0]}%)",
		] )+"<extra></extra>")
	taxp_fig_B.update_traces(
		hovertemplate="<br>".join([
			"<b>Species</b>: %{customdata[2]}",
			"<b>Reads</b>: %{customdata[1]} (%{customdata[0]}%)",
		] )+"<extra></extra>", selector={"name":"<b>Other (grouped)</b>"})

	taxp_fig_B.update_yaxes(showticklabels=False, title_text =None)
	taxp_fig_B.update_xaxes(title_text = '<b> Reads (count) </b>')
	taxp_fig_B.update_layout(yaxis=dict(range=[0.6,1]), 
						plot_bgcolor='white', 
						title_text='<b>'+output+'</b>', 
						title_font=dict(size=20))

	taxp_fig_B_string = pio.to_html(taxp_fig_B, full_html=False)
	taxp_fig_B_data = {
		'taxp_fig_B_'+tool : taxp_fig_B_string,
	}

	return(taxp_fig_B_data)


#
def process_bracken(infile_bracken):
	df_taxp_kraken_bracken = pd.read_csv(infile_bracken, delimiter='\t')

	df_taxp_kraken_bracken = df_taxp_kraken_bracken.sort_values(by='fraction_total_reads', ascending=False)
	df_taxp_kraken_bracken_top_20 = df_taxp_kraken_bracken.nlargest(50,'fraction_total_reads')
	df_taxp_kraken_bracken_other = df_taxp_kraken_bracken[df_taxp_kraken_bracken['name'] != 'unclassified']
	df_taxp_kraken_bracken_other = df_taxp_kraken_bracken[~df_taxp_kraken_bracken.index.isin(df_taxp_kraken_bracken_top_20.index)]
	df_taxp_kraken_bracken_other_frac_total_reads = df_taxp_kraken_bracken_other['fraction_total_reads'].sum()
	df_taxp_kraken_bracken_other_total_reads = df_taxp_kraken_bracken_other['new_est_reads'].sum()
	df_taxp_kraken_bracken_other = pd.DataFrame({'name': ['Other (grouped)'],
						'fraction_total_reads': [df_taxp_kraken_bracken_other_frac_total_reads], 
						'name': ['Other (grouped)'], 'new_est_reads': [df_taxp_kraken_bracken_other_total_reads]})

	df_taxp_kraken_bracken_unclassified = df_taxp_kraken_bracken[df_taxp_kraken_bracken['name'] == 'unclassified']

	df_taxp_kraken_bracken_stacked = pd.concat([df_taxp_kraken_bracken_top_20, df_taxp_kraken_bracken_other, df_taxp_kraken_bracken_unclassified])
	df_taxp_kraken_bracken_stacked['perc_reads'] = (df_taxp_kraken_bracken_stacked['fraction_total_reads']*100).round(2)
	df_taxp_kraken_bracken_stacked['x_value'] = 1

	return(df_taxp_kraken_bracken_stacked)


#
def plot_bracken(infile, tool, output):
	fig_bracken = px.bar(infile,
				x='new_est_reads',
				y='x_value',
				color='name',
				color_discrete_sequence=plot_B_colors,
				orientation='h',
				custom_data=['perc_reads', 'new_est_reads', 'name'],
				labels={'new_est_reads': 'Reads' , 'name': 'Species', 'perc': 'perc_reads'})

	fig_bracken.update_layout(legend_title_text='<b>Species (50 most abundant)</b>')
	fig_bracken.for_each_trace(lambda t: t.update(name = '<i><b>' + t.name +'</i></b>'))
	for t in fig_bracken.data:
		fig_bracken.update_traces(marker={"color": np.where(t.y>0, "lightgrey", t.marker.color)},name = '<b>Other (grouped)</b>', 
							selector={"name":"<i><b>Other (grouped)</i></b>"})

	fig_bracken.update_traces(
		hovertemplate="<br>".join([
			"<b>Species</b>: <i>%{customdata[2]}</i>",
			"<b>Reads</b>: %{customdata[1]} (%{customdata[0]}%)",
		] )+"<extra></extra>"
	)
	fig_bracken.update_traces(
		hovertemplate="<br>".join([
			"<b>Species</b>: %{customdata[2]}",
			"<b>Reads</b>: %{customdata[1]} (%{customdata[0]}%)",
		] )+"<extra></extra>", selector={"name":"<b>Other (grouped)</b>"}
	)

	fig_bracken.update_xaxes(title_text = '<b> Reads (count) </b>')
	fig_bracken.update_yaxes(showticklabels=False, title_text =None)
	fig_bracken.update_layout(yaxis=dict(range=[0.6,1]),
							plot_bgcolor='white',
							title_text='<b>'+output+'</b>',
							title_font=dict(size=20))

	fig_bracken_string = pio.to_html(fig_bracken, full_html=False)
	fig_bracken_data = {
		'fig_bracken_'+tool : fig_bracken_string,
	}

	return(fig_bracken_data)


#
def process_sylph(sylph_infile, sylph_metadata):
	syl_fn_df = pd.read_csv(sylph_metadata, compression='gzip', sep='\t', header=None)
	syl_fn_df.columns = ['Genome_file_1','tax']

	syl = pd.read_csv(sylph_infile, delimiter='\t')
	syl['Genome_file_1'] = syl['Genome_file'].str.split('/').str[-1].str.split("_").str[0:2].apply('_'.join)

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
	syl_merge['rank']
	syl_merge['Species2'] = syl_merge['Species1'].str.replace('s__', '').str.replace('g__', '').str.replace('f__', '').str.replace('o__', '').str.replace('c__', '').str.replace('p__', '').str.replace('d__', '')
	sylph_filtered_df = syl_merge[syl_merge['rank'] == "Species"]
	rcount = sylph_filtered_df['Taxonomic_abundance'].sum()
	top_20 = sylph_filtered_df.nlargest(50,'Taxonomic_abundance')
	other = sylph_filtered_df[sylph_filtered_df['Species2'] != 'unclassified']
	other = sylph_filtered_df[~sylph_filtered_df.index.isin(top_20.index)]
	other_total_reads = other['Taxonomic_abundance'].sum()
	other = pd.DataFrame({'Species2': ['Other (grouped)'], 'Taxonomic_abundance': [other_total_reads]})
	sylph_stacked = pd.concat([top_20, other])
	sylph_stacked['perc'] = ((sylph_stacked['Taxonomic_abundance']/rcount)*100).round(2)
	sylph_stacked['x_value'] = 1

	return(sylph_stacked)


#
def plot_sylph(sylph_stacked, output):
	fig_syl = px.bar(sylph_stacked, 
				x='perc', 
				y='x_value', 
				color='Species2',
				color_discrete_sequence=plot_B_colors,
				orientation='h',
				custom_data=['Species2','perc','Sequence_abundance','Naive_ANI' ],
				labels={'Species2': 'Species2', 'perc': 'perc', 'Sequence_abundance': 'Sequence_abundance','Naive_ANI':'Naive_ANI' })

	fig_syl.update_layout(legend_title_text='<b>Species (50 most abundant)</b>')
	fig_syl.for_each_trace(lambda t: t.update(name = '<i><b>' + t.name +'</i></b>'))

	for t in fig_syl.data:
		fig_syl.update_traces(marker={"color": np.where(t.y>0, "lightgrey", t.marker.color)},
								name = '<b>Other (grouped)</b>', 
								selector={"name":"<i><b>Other (grouped)</i></b>"})

	fig_syl.update_traces(
		hovertemplate="<br>".join([
			"<b>Species</b>: <i>%{customdata[0]}</i>",
			"<b>Taxonomic abundance</b>: %{customdata[1]}%",
			"<b>Sequence abundance</b>: %{customdata[2]}%",
			"<b>Native average nucleotide identity</b>: %{customdata[3]}%",
		] )+"<extra></extra>"
	)
	fig_syl.update_traces(
		hovertemplate="<br>".join([
			"<b>Taxonomic abundance</b>: %{customdata[1]}%",
		] )+"<extra></extra>", selector={"name":"<b>Other (grouped)</b>"}
	)

	fig_syl.update_xaxes(title_text = '<b> Reads (% of assigned) </b>')
	fig_syl.update_yaxes(showticklabels=False, title_text =None)
	fig_syl.update_layout(yaxis=dict(range=[0.6,1]),
						plot_bgcolor='white',
						title_text='<b>'+output+'</b>',
						title_font=dict(size=20))

	fig_syl_string = pio.to_html(fig_syl, full_html=False)
	fig_syl_data = {
		'fig_syl' : fig_syl_string,
	}

	return(fig_syl_data)


#
def merge_ccvals(sample_data, taxp_fig_A_data, taxp_fig_B_data, taxp_fig_A_data_centrifuger, taxp_fig_B_data_centrifuger, fig_bracken_data_kraken2, fig_bracken_data_centrifuger, fig_syl_data):
	legend_count = 0
	context = {}
	context.update(sample_data)

	if taxp_fig_A_data != "None":
		context.update(taxp_fig_A_data)
	if taxp_fig_B_data != "None":
		context.update(taxp_fig_B_data)
		legend_count += 1
		lcount = {
			'kraken2_legend_count': legend_count
		}
		context.update(lcount)
	if fig_bracken_data_kraken2 != "None":
		context.update(fig_bracken_data_kraken2)

	if taxp_fig_A_data_centrifuger != "None":
		context.update(taxp_fig_A_data_centrifuger)
	if taxp_fig_B_data_centrifuger != "None":
		context.update(taxp_fig_B_data_centrifuger)
		legend_count += 1
		lcount = {
			'centrifuger_legend_count': legend_count
		}
		context.update(lcount)
	if fig_bracken_data_centrifuger != "None":
		context.update(fig_bracken_data_centrifuger)

	if fig_syl_data != "None":
		context.update(fig_syl_data)
		legend_count += 1
		lcount = {
			'sylph_legend_count': legend_count
		}
		context.update(lcount)

	return context


#
def render_template(context, template, output):
	template_dir = os.getcwd()
	print(template_dir)
	env = Environment(loader=FileSystemLoader(template_dir))
	template = env.get_template(template)

	html_output = template.render(context)

	html_output = html_output.replace("window.PlotlyConfig = {MathJaxConfig: 'local'};" , "")

	with open(output + '.taxonomy_report.html', 'w') as f:
		f.write(html_output)


# Parse arguments from the command line.
def parse_args():
	description = 'Plots individual and (optionally) merged input from read-based taxonomic classifiers. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--taxpasta_kraken2', required=False, type=str, help='Path to Kraken2 taxpasta output')
	parser.add_argument('--taxpasta_centrifuger', required=False, type=str, help='Path to Centrifuger taxpasta output')
	parser.add_argument('--sylph', type=str, required=False, help='Path to Sylph output')
	parser.add_argument('--bracken_kraken2', required=False, type=str, help='Path to Kraken2 bracken output')
	parser.add_argument('--bracken_centrifuger', required=False, type=str, help='Path to Centrifuger bracken output')
	parser.add_argument('--syl_fn', type=str, required=False, help='Path to gtdb_r214_metadata.tsv.gz file')

	parser.add_argument('--sample_id', required=True, help="Sample ID")
	parser.add_argument('--run_id', required=True, help="Run ID")
	parser.add_argument('--barcode', required=True, help="Barcode")
	parser.add_argument('--sample_type', required=True, help="Sample type")

	parser.add_argument('--logo', required=False, help="Logo")
	parser.add_argument('--report_template', required=True, help="HTML template")

	parser.add_argument('--output', type=str, help='Output file name')
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


###### Main
def main():
	args = parse_args()

	sample_data = process_metadata(args.sample_id, args.run_id, args.barcode, args.sample_type, args.logo)

	if args.taxpasta_kraken2:
		df_taxp, df_taxp_2x1 = process_taxpasta_A(args.taxpasta_kraken2)
		taxp_fig_A_data = taxp_kraken2_plot_A = plot_taxpasta_A(df_taxp_2x1, "Kraken2", args.output)
		df_taxp_stacked = process_taxpasta_B(df_taxp)
		taxp_fig_B_data = taxp_kraken2_plot_B = plot_taxpasta_B(df_taxp_stacked, "Kraken2", args.output)
	else:
		taxp_fig_A_data = "None"
		taxp_fig_B_data = "None"

	if args.taxpasta_centrifuger:
		df_taxp, df_taxp_2x1 = process_taxpasta_A(args.taxpasta_centrifuger)
		taxp_fig_A_data_centrifuger = plot_taxpasta_A(df_taxp_2x1, "Centrifuger", args.output)
		df_taxp_stacked = process_taxpasta_B(df_taxp)
		taxp_fig_B_data_centrifuger = plot_taxpasta_B(df_taxp_stacked, "Centrifuger", args.output)
	else:
		taxp_fig_A_data_centrifuger = "None"
		taxp_fig_B_data_centrifuger = "None"

	if args.bracken_kraken2:
		df_taxp_kraken2_bracken = process_bracken(args.bracken_kraken2)
		fig_bracken_data_kraken2 = plot_bracken(df_taxp_kraken2_bracken, "Kraken2", args.output)
	else:
		fig_bracken_data_kraken2 = "None"

	if args.bracken_centrifuger:
		df_taxp_centrifuger_bracken = process_bracken(args.bracken_centrifuger)
		fig_bracken_data_centrifuger = plot_bracken(df_taxp_centrifuger_bracken, "Centrifuger", args.output)
	else:
		fig_bracken_data_centrifuger = "None"

	if args.sylph and args.syl_fn:
		sylph_results = process_sylph(args.sylph, args.syl_fn)
		fig_syl_data = plot_sylph(sylph_results, args.output)
	else:
		fig_syl_data = "None"

	context = merge_ccvals(sample_data, taxp_fig_A_data, taxp_fig_B_data, taxp_fig_A_data_centrifuger, taxp_fig_B_data_centrifuger, fig_bracken_data_kraken2, fig_bracken_data_centrifuger, fig_syl_data)
	render_template(context, args.report_template, args.output)

if __name__ == "__main__":
	main()
