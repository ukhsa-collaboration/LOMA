#!/usr/bin/env python

# 

__version__ = '0.1'
__date__ = '03-03-2024'
__author__ = 'D.J.BERGER'


###### Imports
import argparse
import os
import tarfile
import random
import string
import shutil
from urllib.request import urlretrieve
import gzip
import re
from git import Repo
import subprocess


###### Functions

def backup_config(config_file):
	backup_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=12))
	backup_file = f"{config_file}.{backup_suffix}"
	shutil.copyfile(config_file, backup_file)


def check_dbdir(database_directory):
	if not database_directory.endswith('/'):
		database_directory += '/'
	if not os.path.exists(database_directory):
		os.makedirs(database_directory)
	if not os.path.isabs(database_directory):
		raise ValueError(f"\n\n\nInput error: The provided path '{database_directory}' is not a full (absolute) path.")

	return database_directory


def rewrite_config(config_file, newfile, parameter):
	search_string = parameter
	with open(config_file, "r") as file:
		lines = file.readlines()

	with open(config_file, "w") as file:
		for line in lines:
			if search_string in line:
				new_line = f'  {search_string} = "{newfile}"\n'
				file.write(new_line)
			else:
				file.write(line)


def get_host_kraken2(database_directory, host_kraken2_url):
	filename = os.path.basename(host_kraken2_url)
	host_kraken2db = os.path.join(database_directory, filename)

	urlretrieve(host_kraken2_url, host_kraken2db)
	with tarfile.open(host_kraken2db, "r:gz") as tar:
		tar.extractall()

	host_kraken2db_path = os.path.join(database_directory, "host_kraken2_db/")

	os.rename("db" , host_kraken2db_path)

	return(host_kraken2db_path)


def get_kraken2db(database_directory, kraken2_url):
	filename = os.path.basename(kraken2_url)
	kraken2db = os.path.join(database_directory, filename)
	kraken2db_dir = re.sub(r'\.tar.gz$', "", kraken2db)

	if not os.path.exists(kraken2db_dir):
		os.makedirs(kraken2db_dir)

	kraken2db_dir += '/'
	urlretrieve(kraken2_url, kraken2db)
	with tarfile.open(kraken2db, "r:gz") as tar:
		tar.extractall(path=kraken2db_dir)

	os.remove(kraken2db)

	return(kraken2db_dir)


def get_centrifugerdb(database_directory, centrifugerdb_url):
	filename = os.path.basename(centrifugerdb_url)
	centrifugerdb_path = os.path.join(database_directory, "centrifuger_db/")
	centrifugerdb_path_full = os.path.join(centrifugerdb_path, filename)

	if not os.path.exists(centrifugerdb_path):
		os.makedirs(centrifugerdb_path)

	for x in range(1, 4):
		centrifugerdb_path_full1 = centrifugerdb_path_full.replace("*", str(x))
		centrifugerdb_url1 = centrifugerdb_url.replace("*", str(x))
		urlretrieve(centrifugerdb_url1, centrifugerdb_path_full1)

	return(centrifugerdb_path)


def get_gtdb_mash(database_directory, gtdb_mash_url):
	filename = os.path.basename(gtdb_mash_url)
	gtdb_mash_path = os.path.join(database_directory, filename)

	urlretrieve(gtdb_mash_url, gtdb_mash_path)

	return(gtdb_mash_path)


def get_checkmdb(database_directory, checkmdb_url):
	filename = os.path.basename(checkmdb_url)
	checkmdb_path = os.path.join(database_directory, filename)
	urlretrieve(checkmdb_url, checkmdb_path)
	checkmdb_dir = re.sub(r'\.tar.gz$', "", checkmdb_path)

	if not os.path.exists(checkmdb_dir):
		os.makedirs(checkmdb_dir)

	checkmdb_dir += '/'

	with tarfile.open(checkmdb_path, "r:gz") as tar:
		tar.extractall(path=checkmdb_dir)

	return checkmdb_dir


def get_sylphdb(database_directory, sylph_url):
	filename = os.path.basename(sylph_url)
	sylphdb_path = os.path.join(database_directory, filename)

	urlretrieve(sylph_url, sylphdb_path)

	return(sylphdb_path)


def get_host_assembly(database_directory, host_assembly_url):
	filename = os.path.basename(host_assembly_url)
	host_assembly_path = os.path.join(database_directory, filename)

	urlretrieve(host_assembly_url, host_assembly_path)

	return(host_assembly_path)


def get_genomad(database_directory, genomad_url):
	filename = os.path.basename(genomad_url)
	genomad_db_path = os.path.join(database_directory, filename)
	urlretrieve(genomad_url, genomad_db_path)
	genomad_dbdir_path = os.path.join(database_directory, "genomad_db/")
	with tarfile.open(genomad_db_path, "r:gz") as tar:
		tar.extractall(path=genomad_dbdir_path)

	os.remove(genomad_db_path)

	return(genomad_dbdir_path)


def get_fdb(database_directory, fdb_url):
	filename = os.path.basename(re.sub(r'\/$', '', fdb_url))
	fdb_path = os.path.join(database_directory, filename)
	fdb_path += '/'

	if not os.path.exists(fdb_path):
		Repo.clone_from(fdb_url, fdb_path)

	return(fdb_path)

def get_gtdb(database_directory, gtdb_url):
	filename = os.path.basename(gtdb_url)
	gtdb_db_path = os.path.join(database_directory, filename)
	urlretrieve(gtdb_url, gtdb_db_path)

	gtdb_dbdir_path = os.path.join(database_directory, "gtdb_db/")

	with tarfile.open(gtdb_db_path, "r:gz") as tar:
		tar.extractall(path=gtdb_dbdir_path)

	os.remove(gtdb_db_path)

	return(gtdb_dbdir_path)


def unpack_repo(database_directory, repodir):
	if not os.path.isfile("virulencefinder_2.0.4--hdfd78af_0.sif"):
		subprocess.run('singularity pull docker://quay.io/biocontainers/virulencefinder:2.0.4--hdfd78af_0', shell=True)

	cwd = os.getcwd()
	os.chdir(repodir)
	subprocess.run('singularity exec ../virulencefinder_2.0.4--hdfd78af_0.sif python INSTALL.py /usr/local/bin/kma', shell=True)

	with open('VERSION','r') as file:
		version = " ".join(line.rstrip() for line in file)

	os.chdir(cwd)

	return(version)

# Parse arguments from the command line
def parse_args():
	description = 'Version: %s, Date: %s, Author: %s' % (__version__, __date__, __author__)
	parser = argparse.ArgumentParser(description=description)

	parser.add_argument("--config_file", required=True, help="Config file")
	parser.add_argument("--db_dir", required=True, help="Database directory")
	parser.add_argument("--genomad", required=False, action='store_true', help="Get geNomad database")
	parser.add_argument("--genomad_url", required=False, default="https://zenodo.org/records/10594875/files/genomad_db_v1.7.tar.gz", help="geNomad database url")
	parser.add_argument("--host_assembly", required=False, action='store_true', help="Get host reference assembly")
	parser.add_argument("--host_assembly_url", required=False, default="https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.fa.gz", help="url of host assembly database url")
	parser.add_argument("--host_kraken2db", required=False, action='store_true', help="Get Host Kraken2 database")
	parser.add_argument("--host_kraken2db_url", required=False, default="https://zenodo.org/records/8339732/files/k2_HPRC_20230810.tar.gz", help="url of host Kraken2 database")
	parser.add_argument("--kraken2db", required=False, action='store_true', help="Get Kraken2 database for taxonomic assignment of reads")
	parser.add_argument("--kraken2db_url", required=False, default="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240605.tar.gz", help="url of Kraken2 database for taxonomic assignment of reads")
	parser.add_argument("--sylphdb", required=False, action='store_true', help="Get Sylph database for taxonomic assignment of reads")
	parser.add_argument("--sylphdb_url", required=False,default="https://storage.googleapis.com/sylph-stuff/gtdb-r220-c200-dbv1.syldb",help="url of Sylph database")
	parser.add_argument("--checkmdb", required=False, action='store_true', help="Get CheckM database")
	parser.add_argument("--checkmdb_url", required=False, default="https://zenodo.org/records/7401545/files/checkm_data_2015_01_16.tar.gz",help="url of Sylph database")
	parser.add_argument("--virulencefinderdb", required=False, action='store_true', help="Get VirulenceFinder database")
	parser.add_argument("--virulencefinderdb_url", required=False, default="https://bitbucket.org/genomicepidemiology/virulencefinder_db/", help="url of VirulenceFinder database")
	parser.add_argument("--pointfinderdb", required=False, action='store_true', help="Get PointFinder database")
	parser.add_argument("--pointfinderdb_url", required=False, default="https://bitbucket.org/genomicepidemiology/pointfinder_db/", help="url of PointFinder database")
	parser.add_argument("--resfinderdb", required=False, action='store_true', help="Get ResFinder database")
	parser.add_argument("--resfinderdb_url", required=False, default="https://bitbucket.org/genomicepidemiology/resfinder_db/", help="url of ResFinder database")
	parser.add_argument("--centrifugerdb", required=False, action='store_true', help="Get Centrifuger database")
	parser.add_argument("--centrifugerdb_url", required=False, default="https://zenodo.org/records/10023239/files/cfr_hpv+gbsarscov2.*.cfr", help="url of Centrifuger database")
	parser.add_argument("--gtdb", required=False, action='store_true', help="Get GTDB-Tk reference database")
	parser.add_argument("--gtdb_url", required=False, default="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz", help="url of GTDB-Tk reference database")
	parser.add_argument("--gtdb_mash", required=False, action='store_true', help="Get GTDB-Tk reference database")
	parser.add_argument("--gtdb_mash_url", required=False, default="https://zenodo.org/records/13731176/files/r220.msh", help="Get GTDB-Tk reference mash database")
	parser.add_argument("--version", action="version", version='Version: %s' % (__version__))

	return(parser.parse_args())


###### Main
def main():
	args = parse_args()

	backup_config(args.config_file)
	database_dir = check_dbdir(args.db_dir)

	if args.centrifugerdb:
		centrifugerdb_path = get_centrifugerdb(database_dir, args.centrifugerdb_url)
		rewrite_config(args.config_file, centrifugerdb_path, "TAXONOMIC_PROFILING.centrifugerdb")

	if args.virulencefinderdb:
		vfdb_path = get_fdb(database_dir, args.virulencefinderdb_url)
		version = unpack_repo(database_dir, vfdb_path)
		rewrite_config(args.config_file, vfdb_path, "VIRULENCEFINDER.db ")
		rewrite_config(args.config_file, version, "VIRULENCEFINDER.db_version")

	if args.resfinderdb:
		rfdb_path = get_fdb(database_dir, args.resfinderdb_url)
		version = unpack_repo(database_dir, rfdb_path)
		rewrite_config(args.config_file, rfdb_path, "RESFINDER.db ")
		rewrite_config(args.config_file, version, "RESFINDER.db_version")

	if args.pointfinderdb:
		pfdb_path = get_fdb(database_dir, args.pointfinderdb_url)
		version = unpack_repo(database_dir, pfdb_path)
		rewrite_config(args.config_file, pfdb_path, "POINTFINDER.db ")
		rewrite_config(args.config_file, version, "POINTFINDER.db_version")

	if args.checkmdb:
		checkmdb_dir_path = get_checkmdb(database_dir, args.checkmdb_url)
		rewrite_config(args.config_file, checkmdb_dir_path, "CHECKM_LINEAGEWF.db")

	if args.kraken2db:
		kraken2db_dir_path = get_kraken2db(database_dir, args.kraken2db_url)
		rewrite_config(args.config_file, kraken2db_dir_path, "TAXONOMIC_PROFILING.krakendb")

	if args.host_kraken2db:
		host_kraken2db_path = get_host_kraken2(database_dir, args.host_kraken2db_url)
		rewrite_config(args.config_file, host_kraken2db_path, "READ_DECONTAMINATION.host_krakendb")

	if args.sylphdb:
		sylphdb_path = get_sylphdb(database_dir, args.sylphdb_url)
		rewrite_config(args.config_file, sylphdb_path, "TAXONOMIC_PROFILING.sylphdb")

	if args.genomad:
		genomad_dbdir_path = get_genomad(database_dir, args.genomad_url)
		rewrite_config(args.config_file, genomad_dbdir_path, "GENOMAD_ENDTOEND.db")

	if args.host_assembly:
		host_assembly_path = get_host_assembly(database_dir, args.host_assembly_url)
		rewrite_config(args.config_file, host_assembly_path, "READ_DECONTAMINATION.host_assembly")

	if args.gtdb:
		gtdb_path = get_gtdb(database_dir, args.gtdb_url)
		rewrite_config(args.config_file, gtdb_path, "GTDBTK_CLASSIFYWF.gtdb_db")

	if args.gtdb_mash:
		gtdb_mash_path = get_gtdb_mash(database_dir, args.gtdb_mash_url)
		rewrite_config(args.config_file, gtdb_mash_path, "GTDBTK_CLASSIFYWF.mash_db")

if __name__ == "__main__":
	main()
