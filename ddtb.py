#!/usr/bin/env python

"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.2
CREATED: 27 March 2019
REVISION:
    -01 April 2019 handle arguments, inputs and outputs
    -02 April 2019 test distance matrix calculation and make subtasks modular
TODO
	-...
    - Absolute paths (args.short1 = os.path.abspath(args.short1))

================================================================
END_OF_HEADER
================================================================
"""

import os
import pandas as pd
import argparse
from misc import check_file_exists, import_to_pandas
import ddtb_add
import ddtb_compare
import ddtb_update
import ddtb_extract


def get_arguments():

    #Define parser and program
    parser = argparse.ArgumentParser(prog = "VCF_DDTB.py", description= "VCF_DDTB manages a custom snp database") 

    #Define subtask/subparsers
    subparsers = parser.add_subparsers( dest = "subtask", help = "new / update / compare / extract commands either add new samples, compare or discard exixting samples")

    new_parser = subparsers.add_parser("new", help = "Create new ddbb with presence/absence of snp")
    update_parser = subparsers.add_parser("update", help = "Add new sample using a list of variants, files supplied or files on folder")
    compare_parser = subparsers.add_parser("compare", help = "Comapare samples supplied or all samples to obtain a pirwise matrix")
    extract_parser = subparsers.add_parser("extract", help = "Remove samples supplied from databse")


    #new_parser
    new_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the newd ddbb")
    new_parser.add_argument("-v", "--vcf", required= False, action='store_true', help="Database will use vcf files instead of ")
    new_parser.add_argument("-s", "--suffix", required= False, type=str, default=".SNP.final.vcf", help="Suffix to filter within vcf with similar suffix")
    new_parser.add_argument("-r", "--recalibrate", required= False, type=str, default=False, help="Folder with tab files to asses discrepancies after comparing")


    new_exclusive = new_parser.add_mutually_exclusive_group(required= True)

    new_exclusive.add_argument("-F", "--folder",  dest = "folder", metavar="folder",required= False, type=str, help="Folder containinig files with snp positions")
    new_exclusive.add_argument("-f", "--file",  dest = "single_file", metavar="file[s]", required= False, type=str, help="individual files with snp positions")
    new_exclusive.add_argument("-l", "--list",  dest = "snp_list", metavar="list", required= False, help="file with a list of positions in a column, not vcf format")


    #update_parser
    update_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")
    update_parser.add_argument("-s", "--snp-final", required= False, action='store_true', help="Database will use snp.fila instead of gatk vcf ")
    update_parser.add_argument("-d", "--database",  dest = "update_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")

    update_exclusive = update_parser.add_mutually_exclusive_group()

    update_exclusive.add_argument("-F", "--folder",  dest = "folder", metavar="folder",required= False, type=str, help="Folder containinig files with snp positions")
    update_exclusive.add_argument("-f", "--file",  dest = "single_file", metavar="file[s]", required= False, type=str, help="individual files with snp positions")
    update_exclusive.add_argument("-l", "--list",  dest = "snp_list", metavar="list", required= False, help="file with a list of positions in a column, not vcf format")

    update_parser.add_argument("-b", "--backup",  dest = "backup", action="store_true", help="Creates an aditional database with the date as backup")

    #compare_parser
    compare_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
    compare_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= False, metavar="filename", help="REQUIRED: file name, including PATH, of the matrix comparison")

    compare_exclusive = compare_parser.add_mutually_exclusive_group()

    compare_exclusive.add_argument("-a", "--all",  dest = "all_compare", action="store_true", required= False, help="All files in supplied database will be compared")
    compare_exclusive.add_argument("-s", "--samples",  dest = "samples_compare", metavar="sample_name[s]", nargs="+", required= False, help="Sample names supplied will be compared")

    #extract_parser
    extract_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database with the sample to remove")
    extract_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")

    extract_parser.add_argument("-s", "--samples",  dest = "samples_extract", metavar="sample name[s]", nargs="+", required= True, help="Sample names supplied will be removed")

    parser.add_argument("--version", action="version", version="%(prog)s 0.1")

    arguments = parser.parse_args()

    return arguments


 #TODO: 
    #- better implementation of sample name
    #- exit if file don't exist
    #- create a file with the date as backup
    #- handle when sample exist
    #- handle cwd instead of full path
    #- if __name__ == "__main__":



args = get_arguments()

print(args)

if args.subtask == 'new' :
    ddtb_add.ddtb_add(args)
elif args.subtask == 'compare':
    ddtb_compare.ddtb_compare(args)
elif args.subtask == 'update':
    ddtb_add.ddtb_add(args)
    #ddtb_update.ddtb_update(args)
elif args.subtask == 'extract':
    ddtb_extract.ddtb_extract(args)



