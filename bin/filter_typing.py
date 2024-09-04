#!/usr/bin/env python

# Merges results from multiple E. coli/Shigella typing tools to provide a summarized output and a estimation of assignment/species/pathotype.
# Note: Hard-coded cluter IDs in places due to the way various tools function. These may need to be updated at some point.

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports
import pandas as pd
import argparse

###### Functions


def process_mlst(mlst_results, krocus_results):
        mlst_designation = []
        mlst = pd.read_csv(mlst_results, sep='\t')
        mlst = mlst[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
        mlst.columns = ['Name', "MLST: Sequence type", "Internal UKHSA: Clonal complex (mlst)"]
        krocus = pd.read_csv(krocus_results, sep='\t')
        krocus = krocus[['Name', "Sequence type", "Internal UKHSA: Clonal complex"]]
        krocus.columns = ['Name', "Krocus: Sequence type", "Internal UKHSA: Clonal complex (Krocus)"]

        cc_merged = pd.merge(mlst, krocus,on=['Name'], how = 'outer')

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

        mlst_designation = cc_merged['Internal UKHSA: Clonal complex'][0]

        return(mlst_designation, cc_merged)

def process_salmonella(seqsero2_results, sistr_results):
        seqsero2 = pd.read_csv(seqsero2_results, sep='\t')
        sistr = pd.read_csv(sistr_results, sep='\t')

        if "Paratyphi B" in seqsero2.at[0, 'Predicted serotype']:
                salmonella_designation = "Paratyphi_B"
        else:
                salmonella_designation = "Typhi"

        return(salmonella_designation, seqsero2)


# Write an output file with tpying info as part of file name
def save_results(designation, file, bin, output):
        file.to_csv(str(output) + "." + str(bin) + "." + str(designation) + ".ftype.tsv", encoding='utf-8', index=False, sep='\t', na_rep='NA')


# Parse arguments from the command line.
def parse_args():
        description = 'Merge MLST and Krocus results, get a consensus sequence type (ST) and clonal complex (CC), output results file with ST appended. Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
        parser = argparse.ArgumentParser(description="description")
        parser.add_argument('--mlst', required=False, help="mlst results")
        parser.add_argument('--krocus', required=False, help="Krocus results")
        parser.add_argument('--seqsero2', required=False, help="SeqSero2 results")
        parser.add_argument('--sistr', required=False, help="SISTR results")
        parser.add_argument('--output', required=True, help="Output TSV name")
        parser.add_argument('--bin', required=True, help="Bin name")

        parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

        return(parser.parse_args())


###### Main
def main():
        args = parse_args()

        if args.mlst and args.krocus:
                mlst_designation, cc_merged = process_mlst(args.mlst, args.krocus)
                save_results(mlst_designation, cc_merged, args.bin, args.output)

        elif args.seqsero2 and args.sistr:
                salmonella_designation, seqsero2 = process_salmonella(args.seqsero2, args.sistr)
                save_results(salmonella_designation, seqsero2, args.bin, args.output)

        else:
                print("Incorrect input, should be either: --seqsero2 & --sistr or --mlst & --krocus")
                exit

if __name__ == "__main__":
        main()
