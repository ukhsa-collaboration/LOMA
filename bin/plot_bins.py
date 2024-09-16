#!/usr/bin/env python

# Plots a graphical summary of metagenome assembled genomes in interactive HTML format with some additional assembly quality control plots.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports

import math
import base64
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import plotly.express as px
import argparse
from jinja2 import Environment, FileSystemLoader
from io import BytesIO
import plotly.io as pio
import os

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


def process_bintax(bin_tax):
	bin_tax= pd.read_csv(bin_tax, sep='\t')
	bin_tax['Bin'] = bin_tax['user_genome'].str.split('.').str[1]
	bin_tax['Species1'] = bin_tax['classification'].str.replace(
		r'd__;p__;o__;f__;g__;s__$', '', regex=True).str.replace(
			r';p__;o__;f__;g__;s__$', '', regex=True).str.replace(
				r';c__;o__;f__;g__;s__$', '', regex=True).str.replace(
					r';o__;f__;g__;s__$', '', regex=True).str.replace(
						r';f__;g__;s__$', '', regex=True).str.replace(
							r';g__;s__$', '', regex=True).str.replace(
								r';s__$', '', regex=True).str.split(';').str[-1]

	bin_tax['rank'] = np.select(
		[bin_tax['Species1'].str.contains('s__'),
		bin_tax['Species1'].str.contains('g__'),
		bin_tax['Species1'].str.contains('f__'),
		bin_tax['Species1'].str.contains('o__'),
		bin_tax['Species1'].str.contains('c__'),
		bin_tax['Species1'].str.contains('p__'),
		bin_tax['Species1'].str.contains('d__'),
		bin_tax['Species1'].str.contains('Unclassified')],
		['Species', 'Genus','Family','Order','Class','Phylum','Domain','Unclassified'],
		default=np.nan,
	)

	bin_tax['Species2'] = bin_tax['Species1'].str.replace(
		's__', '').str.replace(
			'g__', '').str.replace(
				'f__', '').str.replace(
					'o__', '').str.replace(
						'c__', '').str.replace(
							'p__', '').str.replace(
								'd__', '')
	return(bin_tax)


def process_skani(skani, gtdb_fn):
	skani_contigs = pd.read_csv(skani, delimiter='\t')
	gtdb_fn_df = pd.read_csv(gtdb_fn, compression='gzip', sep='\t', header=None)
	gtdb_fn_df.columns = ['Genome_file_1','tax']

	skani_contigs['Genome_file_1'] = skani_contigs['Ref_file'].str.split('/').str[-1].str.split("_").str[0:2].apply('_'.join)
	skani_contigs_merge = pd.merge(skani_contigs, gtdb_fn_df,on=['Genome_file_1'], how = 'left')

	skani_contigs_merge['Species1'] = skani_contigs_merge['tax'].str.replace(
		r'd__;p__;o__;f__;g__;s__$', '', regex=True).str.replace(
			r';p__;o__;f__;g__;s__$', '', regex=True).str.replace(
				r';c__;o__;f__;g__;s__$', '', regex=True).str.replace(
					r';o__;f__;g__;s__$', '', regex=True).str.replace(
						r';f__;g__;s__$', '', regex=True).str.replace(
							r';g__;s__$', '', regex=True).str.replace(
								r';s__$', '', regex=True).str.split(';').str[-1]


	skani_contigs_merge['rank'] = np.select(
		[skani_contigs_merge['Species1'].str.contains('s__'),
		skani_contigs_merge['Species1'].str.contains('g__'),
		skani_contigs_merge['Species1'].str.contains('f__'),
		skani_contigs_merge['Species1'].str.contains('o__'),
		skani_contigs_merge['Species1'].str.contains('c__'),
		skani_contigs_merge['Species1'].str.contains('p__'),
		skani_contigs_merge['Species1'].str.contains('d__'),
		skani_contigs_merge['Species1'].str.contains('Unclassified')],
		['Species', 'Genus','Family','Order','Class','Phylum','Domain','Unclassified'], default=np.nan)

	skani_contigs_merge['Species2'] = skani_contigs_merge['Species1'].str.replace(
		's__', '').str.replace(
			'g__', '').str.replace(
				'f__', '').str.replace(
					'o__', '').str.replace(
						'c__', '').str.replace(
							'p__', '').str.replace('d__', '')

	skani_contigs_merge2 = skani_contigs_merge.groupby("Query_name").first()
	skani_contigs_merge2 = skani_contigs_merge2.reset_index()
	skani_contigs_merge3 = skani_contigs_merge2[['Query_name', 'Species2','rank', 'ANI','Align_fraction_ref','Align_fraction_query','Ref_name']]
	skani_contigs_merge3.columns = ['name_y','contig_species2','contig_rank','contig_ANI','contig_align_fraction_ref','contig_align_fraction_query','Reference match']

	return(skani_contigs_merge3)


def process_checkm(checkm_metrics):
	checkm = pd.read_csv(checkm_metrics, sep='\t')
	checkm.columns = ['name','Marker_lineage','genomes','markers','marker_sets','0','1','2','3','4','5+','Completeness','Contamination','Strain heterogeneity']
	checkm['Bin'] = checkm['name'].str.split('.').str[1]

	conditions = [
		(checkm['Completeness'] == 'NaN'),
		(checkm['Completeness'] >= 90) & (checkm['Contamination'] <= 5),
		(checkm['Completeness'] >= 90) & (checkm['Contamination'] > 5) & (checkm['Contamination'] <= 10),
		(checkm['Completeness'] >= 70) & (checkm['Completeness'] < 90) & (checkm['Contamination'] <= 10),
		(checkm['Completeness'] >= 50) & (checkm['Completeness'] < 70) & (checkm['Contamination'] <= 10),
		(checkm['Completeness'] <= 50),
		(checkm['Contamination'] > 10)]

	values = ['NA', 'High quality', 'Medium quality', 'Medium quality', 'Partial assembly', 'QC fail', 'QC fail']
	checkm['bin_qual'] = np.select(conditions, values)
	checkm['bin_qual'] = checkm['bin_qual'].str.replace('0','TBD')

	return(checkm)


def merge_stats(asm_metrics, fstat_metrics, cov_metrics, plasmid_metrics, bintax_metrics, skani_metrics, checkm_metrics):
	asm_stats = pd.read_csv(asm_metrics)
	asm_stats.columns = ['sampleid','Bin','assembly_length_bp','scaffold_count','scaffold_N50_bp','scaffold_N90_bp','contig_count','contig_N50_bp','contig_N90_bp','GC_perc','gaps_count','gaps_sum_bp','gaps_perc']

	fstats = pd.read_csv(fstat_metrics, sep='\t', header=None)
	fstats.columns = ['name','len','GC','N_count']

	cov = pd.read_csv(cov_metrics, sep='\t')
	cov.columns = ['name','startpos','endpos','numreads','covbases','coverage','meandepth','meanbaseq','meanmapq']

	plasmid_score = pd.read_csv(plasmid_metrics, delimiter='\t')
	plasmid_score.columns = ['name', 'length','topology', 'n_genes', 'genetic_code', 'plasmid_score','fdr','n_hallmarks', 'marker_enrichment', 'conjugation_genes', 'amr_genes']
	plasmid_score2 = plasmid_score[["name", "plasmid_score"]]

	m1x = pd.merge(fstats, plasmid_score2,on=['name'], how = 'outer' )
	m1 = pd.merge(m1x, cov,on=['name'])
	m1['Bin'] = m1['name'].str.split('.').str[1]
	m2 = pd.merge(m1, bintax_metrics,on=['Bin'], how = 'outer')

	m3 = pd.merge(checkm_metrics, m2,on=['Bin'], how = 'outer')
	m3b = pd.merge(asm_stats, m3, on=['Bin'], how = 'outer')

	m3b['Bin2'] = m3b['Bin'].replace("bin_","", regex=True).replace("unbinned","Unbinned", regex=True)
	m3b['rGC'] = m3b['GC'].round(1)
	m3b['rmeandepth'] = m3b['meandepth'].round(1)
	m3b['rcoverage'] = m3b['coverage'].round(1)
	m3b['rlen'] = (m3b['len']/1000).round(1)
	m3b['Bin3'] = m3b['Bin2'].astype(str) + (" (") + m3b['Species2'] + (")")
	m3b.Bin3 = m3b.Bin3.fillna('Unbinned')

	m6 = m3b.fillna('N/A')
	m7 = pd.merge(m6, skani_metrics, on=['name_y'], how = 'outer')
	m7['contig_species3'] = m7['contig_species2'].fillna('N/A')
	m7['contig_rank2'] = m7['contig_rank'].fillna('N/A')
	m7['contig_ANI2'] = m7['contig_ANI'].fillna('N/A')
	m8 = m7.sort_values(by='Bin')

	return(m7, m8, asm_stats)


def plot_bins(mq_out, output):
	bin_dist = px.scatter(mq_out,
			x='GC',
			y='meandepth',
			opacity=.75,
			size='len',
			log_y=True,
			size_max=40,
			color="Bin2",
			color_discrete_sequence=["#ff4a24","#19D3F3","#F58518","#e60019","#006fdf","#FECB52","#636EFA","#80ba5a","#b596ed","#009b3b","#bc80bd","#882255","#117733","#88ccee","#ff0b6e","#00d395","#f2b701","#999933","#E45756","#332288","#cf1c90","#FFA15A","#ccebc5","#80b1d3","#FF9DA6","#002b90","#008695","#11a579","#3969ac","#44aa99","#54A24B","#72B7B2","#7f3c8d","#9D755D","#aa4499","#AB63FA","#B279A2","#b3de69","#B6E880","#bebada","#c44a86","#cc6677","#ddcc77","#de57ff","#e59800","#e68310","#e73f74","#EECA3B","#EF553B","#fb8072","#fccde5","#FF6692","#FF97FF"],
			custom_data=['name_y', 'rlen', 'rGC','rcoverage','rmeandepth',
							'Bin2','Species2','rank','bin_qual','closest_genome_ani',
							'Completeness','Contamination','Strain heterogeneity','contig_species3',
						'contig_rank2','contig_ANI2','contig_align_fraction_ref','contig_align_fraction_query', 'plasmid_score'],
			)

	for t in bin_dist.data:
		bin_dist.update_traces(marker={"color": ("lightgrey")},selector={"name":"Unbinned"})

	bin_dist.for_each_trace(lambda t: t.update(marker_line_color = t.marker.color, opacity=1))

	bin_dist.for_each_trace(lambda t: t.update(name = '<b>' + t.name +'</b>' + " (" + '<i>' + f"{t.customdata[0][6]}" + '</i>' + ")"))

	bin_dist.for_each_trace(lambda t: t.update(name = '<b>' + "Unbinned" +'</b>'), 	selector=dict(name = "<b>Unbinned</b> (<i>N/A</i>)"))

	bin_dist.update_traces(
		hovertemplate="<br>".join([
			"<b>Contig name</b>: %{customdata[0]}",
			"<b>Length </b>: %{customdata[1]} kbp",
			"<b>GC</b>: %{customdata[2]}%",
			"<b>Mean depth</b>: %{customdata[4]}",
			"<b>Contig coverage</b>: %{customdata[3]}%",
			"<b>Contig taxonomy (rank)</b>: <i>%{customdata[13]}</i> (%{customdata[14]})",
			"<b>Contig average nucleotide identity</b>: %{customdata[15]}%",
			"<b>Plasmid score</b>: %{customdata[18]}",
			"<b>--</b>",
			"<b>Assigned bin</b>: %{customdata[5]}",
			"<b>Assigned taxonomy (rank)</b>: <i>%{customdata[6]}</i> (%{customdata[7]})",
			"<b>Average nucleotide identity</b>: %{customdata[9]}%",
			"<b>Completeness</b>: %{customdata[10]}% (%{customdata[8]})",
			"<b>Contamination (Strain heterogenity)</b>: %{customdata[11]}% (%{customdata[12]})",
		] )+"<extra></extra>"
	)

	bin_dist.update_traces(selector={"name":"<b>Unbinned</b>"},
		hovertemplate="<br>".join([
			"<b>Contig name</b>: %{customdata[0]}",
			"<b>Length </b>: %{customdata[1]} kbp",
			"<b>GC</b>: %{customdata[2]}%",
			"<b>Mean depth</b>: %{customdata[4]}",
			"<b>Contig coverage</b>: %{customdata[3]}%",
			"<b>Contig taxonomy (rank)</b>: <i>%{customdata[13]}</i> (%{customdata[14]})",
			"<b>Contig average nucleotide identity</b>: %{customdata[15]}%",
			"<b>Plasmid score</b>: %{customdata[18]}",
			"<b>--</b>",
			"<b>Assigned bin</b>: Unbinned",
		] )+"<extra></extra>"
	)

	bin_dist.update_xaxes(title_text = '<b> GC (%) </b>', showline=True, linewidth=1, linecolor='black')
	bin_dist.update_yaxes(title_text ='<b> Log\u2081\u2080(mean coverage) </b>', showline=True, linewidth=1, linecolor='black')
	bin_dist.update_layout(title_text='<b>'+output+'</b>', title_font=dict(size=20), plot_bgcolor='white', legend_title_text='<b>Assigned bin</b>', height=850, width=1500)

	bin_dist_string = pio.to_html(bin_dist, full_html=False)
	bin_dist_data = {
		'bin_dist' : bin_dist_string,
	}

	return(bin_dist_data)


def assembly_summary(m7):
	m8 = m7[m7["Contamination"].astype(str).str.contains("N/A")!=True]
	mincompl = m8['Completeness'].min()
	maxcont = m8['Contamination'].max()
	roundedqual = (math.floor(mincompl/5)*5)
	roundedcont = (math.ceil(maxcont/10)*10)
	maxcontig = m8['contig_count'].max()
	roundedcontig = (math.ceil(maxcontig/20)*20)
	m9 = m8[["Completeness", "Contamination", "bin_qual","contig_count" ]]
	m10 = m9.drop_duplicates()

	m11 = m10['bin_qual'].value_counts().rename_axis('bin_qual').reset_index(name='counts')

	return(m11, m10, roundedqual, roundedcont)


def plot_assembly_summary(m10, roundedqual, roundedcont):
	palette = {'High quality': '#7f58af',
								'Medium quality': "#64c5eb",
								'Partial assembly': "#e83d8a",
								"QC fail": "#feb326"}
	gfew = sns.jointplot(
		data=m10,
		x="Completeness",
		alpha=.85,
		hue="bin_qual",
		palette=palette,
		s=16,
		y="Contamination",
		hue_order = ['High quality', 'Medium quality', 'Partial assembly','QC fail'],
		kind='scatter',
		xlim = (roundedqual-2.5,100), ylim = (-1,roundedcont))

	gfew.set_axis_labels("Completeness (%)","Contamination (%)", fontweight='bold')
	gfew.figure.tight_layout()
	gfew.ax_joint.legend(loc='upper left')

	gfew.savefig('figaaa1.png', bbox_inches="tight", dpi=699)


def plot_assembly_stats(m11, asm_stats):
	gridspec = dict(height_ratios=[0.12,0.01, 0.025, 0.025, 0.025,0.025])
	palette = {'High quality': '#7f58af',
								'Medium quality': "#64c5eb",
								'Partial assembly': "#e83d8a",
								"QC fail": "#feb326"}
	fig, axs = plt.subplots(ncols=1, nrows=6, gridspec_kw=gridspec, figsize = (7, 13))

	qqq = sns.stripplot(x="contig_count", data=asm_stats, jitter = True, color="grey",alpha=.75, ax=axs[2])
	qqq = sns.boxplot(x="contig_count", data=asm_stats, showfliers=False, fill=False, color="black", ax=axs[2])
	qqq.set_xlabel("Number of contigs", fontdict={'weight': 'bold'})
	qqq.set_ylabel("")
	qqq.tick_params(left=False)

	qqq2 = sns.stripplot(x="contig_N50_bp", data=asm_stats, jitter = True, color="grey",alpha=.75, ax=axs[3])
	qqq2 = sns.boxplot(x="contig_N50_bp", data=asm_stats, showfliers=False, fill=False, color="black", ax=axs[3])
	qqq2.set_xlabel(r'Contig N$\bf{_{50}}$ (Mb)', fontdict={'weight': 'bold'})
	qqq2.set_ylabel("")
	qqq2.tick_params(left=False)

	qqq3 = sns.stripplot(x="assembly_length_bp", data=asm_stats, jitter = True, color="grey",alpha=.75, ax=axs[4])
	qqq3 = sns.boxplot(x="assembly_length_bp", data=asm_stats, showfliers=False, fill=False, color="black", ax=axs[4])
	qqq3.set_xlabel("Assembly length (Mb)", fontdict={'weight': 'bold'})
	qqq3.set_ylabel("")
	qqq3.tick_params(left=False)

	qqq4 = sns.stripplot(x="GC_perc", data=asm_stats, jitter = True, color="grey",alpha=.75, ax=axs[5])
	qqq4 = sns.boxplot(x="GC_perc", data=asm_stats, showfliers=False, fill=False, color="black", ax=axs[5])
	qqq4.set_xlabel("GC (%)", fontdict={'weight': 'bold'})
	qqq4.tick_params(left=False)
	sns.despine( ax=axs[5], left=True)
	sns.despine( ax=axs[4], left=True)
	sns.despine( ax=axs[3], left=True)
	sns.despine( ax=axs[2], left=True)
	sns.despine(left=False, ax=axs[0])

	qqq0 = sns.barplot(y="counts",x="bin_qual", data=m11, ax=axs[0], hue="bin_qual", palette=palette)
	qqq0.set_xticklabels(qqq0.get_xticklabels(),rotation=30)
	qqq0.set_xlabel("", fontdict={'weight': 'bold'}, rotation=45  )
	qqq0.set_ylabel("Count", fontdict={'weight': 'bold'})
	qqq0 = sns.barplot(y="counts",x="bin_qual", data=m11, hue="bin_qual", ax=axs[1])
	axs[1].set_visible(False)

	fig.subplots_adjust(hspace=1)
	fig.savefig('figaaa.png', bbox_inches="tight")


def merge_plots(m7, output):

	f, axarr = plt.subplots(1, 2, figsize=(14, 7))

	axarr[0].imshow(mpimg.imread('figaaa1.png'))
	axarr[1].imshow(mpimg.imread('figaaa.png'))

	[ax.set_axis_off() for ax in axarr.ravel()]

	plt.tight_layout()

	buffer = BytesIO()
	plt.savefig(buffer, bbox_inches="tight", format="png")
	buffer.seek(0)
	image_base64 = base64.b64encode(buffer.read()).decode('utf-8')
	buffer.close()

	fig_data = {
		'image_base64' : image_base64,
	}


	f.savefig('figaa2.png', bbox_inches="tight")

	outstringtab = output + '.bin_summary.tsv'
	outstringcontigtab = output + '.contig_summary.tsv'

	outtab = m7[["Bin","GC_perc","assembly_length_bp","contig_count","contig_N50_bp","scaffold_count","scaffold_N50_bp","gaps_count","gaps_sum_bp","Marker_lineage","Completeness","Contamination", "Strain heterogeneity","classification","closest_genome_ani","closest_genome_af","Species2","rank","bin_qual"]].drop_duplicates()
	outtab = outtab[~outtab['Bin'].isin(['unbinned'])]
	outtab = outtab.replace('N/A','')

	outtab3 = outtab
	outtab2 = outtab.rename(columns={'GC_perc': 'GC (%)','assembly_length_bp': 'Assembly length (bp)', "contig_count":"Contigs (count)", "contig_N50_bp":"Contig N50 (bp)", "scaffold_count":"Scaffolds (count)", "scaffold_N50_bp":"Scaffold N50 (bp)", "gaps_count":"Gaps (count)", "gaps_sum_bp":"Total length N's (bp)", "Marker_lineage":"CheckM: Marker_lineage", "Completeness":"CheckM: Completeness", "Contamination":"CheckM: Contamination", "Strain heterogeneity":"CheckM: Strain heterogeneity", "classification":"GTDB-TK: classification", "closest_genome_ani":"GTDB-TK: FastANI average nucleotide identity", "closest_genome_af":"GTDB-TK: FastANI allele frequency", "Species2":"Inferred taxononomic classification", "rank":"Inferred classification rank", "bin_qual":"Bin quality classification"})

	outcont_tab = m7[["name_y","len","GC","N_count","coverage","meandepth","meanmapq","Reference match","contig_align_fraction_ref","contig_align_fraction_query","contig_species3","contig_rank2","contig_ANI2"]]
	outcont_tab = outcont_tab.replace('N/A','')
	outcont_tab2 = outcont_tab
	outcont_tab = outcont_tab.rename(columns={"name_y":"Contig name","len":"Contig length (bp)","GC":"Contig GC (%)","N_count":"Contig N's (count)","coverage":"Contig coverage (%)","meandepth":"Mean depth","meanmapq":"Mean mapping quality","Reference match":"skani: Nearest reference match","contig_align_fraction_ref":"skani: Fraction of reference covered by alignments","contig_align_fraction_query":"skani: Fraction of query covered by alignments","contig_species3":"skani: Inferred taxononomic classification","contig_rank2":"skani: Inferred taxononomic classification rank","contig_ANI2":"skani: Contig average nucleotide identity to reference"})


	outtab2.to_csv(outstringtab, encoding='utf-8', index=False, sep='\t')
	outcont_tab.to_csv(outstringcontigtab, encoding='utf-8', index=False, sep='\t')

	outtab3 = outtab3[["Bin","GC_perc","assembly_length_bp","scaffold_count","scaffold_N50_bp","gaps_count","Completeness","Contamination", "Strain heterogeneity","Species2","rank","bin_qual"]]
	outtab3['assembly_length_bp'] = outtab3['assembly_length_bp'].astype(int)
	outtab3['scaffold_count'] = outtab3['scaffold_count'].astype(int)
	outtab3['gaps_count'] = outtab3['gaps_count'].astype(int)
	outtab3['scaffold_N50_bp'] = outtab3['scaffold_N50_bp'].astype(int)
	outtab3['Bin'] = outtab3['Bin'].str.replace('bin_', '')
	binning_merged_data = outtab3.to_dict(orient='records')

	outcont_tab2 = outcont_tab2[["name_y","len","GC","N_count","coverage","meandepth","meanmapq","Reference match","contig_species3","contig_rank2","contig_ANI2"]]
	outcont_tab2['len'] = outcont_tab2['len'].astype(int)

	outcont_tab2_data = outcont_tab2.to_dict(orient='records')

	return(binning_merged_data, fig_data, outcont_tab2_data)


def merge_ccvals(sample_data, bin_dist_data, binning_merged_data, fig_data, outcont_tab2_data):
	legend_count = 0
	context = {}
	context.update(sample_data)
	context.update(bin_dist_data)
	context['binning_merged_data'] = binning_merged_data
	context['outcont_tab2_data'] = outcont_tab2_data
	context.update(fig_data)

	return context


#
def render_template(context, template, output):
	template_dir = os.getcwd()
	env = Environment(loader=FileSystemLoader(template_dir))
	template = env.get_template(template)

	html_output = template.render(context)

	html_output = html_output.replace("window.PlotlyConfig = {MathJaxConfig: 'local'};" , "")
	html_output = html_output.replace("<em>Unclassified</em>","Unclassified")
	html_output = html_output.replace("<em>Unclassified Bacteria</em>","Unclassified Bacteria")
	html_output = html_output.replace("<em></em> ()","")

	with open(output + '.summary_binning_report.html', 'w') as f:
		f.write(html_output)




def parse_args():
	description = 'Plot aggregate metagenome assembled genome qc and analysis results. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('--asm_stats', type=str, help='Path to assembly stats file', required=True)
	parser.add_argument('--cov', type=str, help='Path to coverage file', required=True)
	parser.add_argument('--fstats', type=str, help='Path to contig stats file', required=True)
	parser.add_argument('--bin_tax', type=str, help='Path to GTDB-tk output  file', required=True)
	parser.add_argument('--checkm', type=str, help='Path to CheckM file', required=True)
	parser.add_argument('--skani', type=str, help='Path to Skani results file', required=True)
	parser.add_argument('--genomad_plasmid', type=str, help='Path to geNomad plasmid summary file', required=True)
	parser.add_argument('--gtdb_fn', type=str, help='Path to gtdb_r214_metadata.tsv.gz file', required=True)

	parser.add_argument('--sample_id', required=True, help="Sample ID")
	parser.add_argument('--run_id', required=True, help="Run ID")
	parser.add_argument('--barcode', required=True, help="Barcode")
	parser.add_argument('--sample_type', required=True, help="Sample type")
	parser.add_argument('--logo', required=False, help="Logo")
	parser.add_argument('--report_template', required=True, help="HTML template")

	parser.add_argument('--output', type=str, help='Output file name', required=True)
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return parser.parse_args()


def main():
	args = parse_args()
	sample_data = process_metadata(args.sample_id, args.run_id, args.barcode, args.sample_type, args.logo)
	bintax_metrics = process_bintax(args.bin_tax)
	skani_metrics = process_skani(args.skani, args.gtdb_fn)
	checkm_metrics = process_checkm(args.checkm)
	m7, mq_out, asm_stats = merge_stats(args.asm_stats, args.fstats, args.cov, args.genomad_plasmid, bintax_metrics, skani_metrics, checkm_metrics)
	bin_dist_data = plot_bins(mq_out, args.output)
	m11, m10, roundedqual, roundedcont = assembly_summary(mq_out)
	plot_assembly_summary(m10, roundedqual, roundedcont)
	plot_assembly_stats(m11, asm_stats)
	binning_merged_data, fig_data, outcont_tab2_data = merge_plots(m7, args.output)

	context = merge_ccvals(sample_data, bin_dist_data, binning_merged_data, fig_data, outcont_tab2_data)
	render_template(context, args.report_template, args.output)


if __name__ == "__main__":
	main()
