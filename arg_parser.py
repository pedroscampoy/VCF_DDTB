import argparse

#Define parser and program
parser = argparse.ArgumentParser(prog = "VCF_DDTB.py", description= "VCF_DDTB manages a custom snp database") 

#Define subtask/subparsers
subparsers = parser.add_subparsers( dest = "subtask", help = "update / compare commands either update or compare exixting samples",)

update_parser = subparsers.add_parser("update", help = "Add new sample using a list of variants, files supplied or files on folder")
compare_parser = subparsers.add_parser("compare", help = "Comapare samples supplied or all samples to obtain a pirwise matrix")
extract_parser = subparsers.add_parser("extract", help = "Remove samples supplied from databse")


#update_parser
update_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
update_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")


update_exclusive = update_parser.add_mutually_exclusive_group()

update_exclusive.add_argument("-F", "--folder",  dest = "folder", metavar="folder",required= False, help="Folder containinig files with snp positions")
update_exclusive.add_argument("-f", "--file",  dest = "single_file", metavar="file[s]", required= False, help="individual files with snp positions")
update_exclusive.add_argument("-l", "--list",  dest = "snp_list", metavar="list", required= False, help="file with a list of positions in a column, not vcf format")

update_parser.add_argument("-b", "--backup",  dest = "backup", action="store_true", help="Creates an aditional database with the date as backup")


#compare_parser
compare_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
compare_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the matrix comparison")


compare_exclusive = compare_parser.add_mutually_exclusive_group()

compare_exclusive.add_argument("-a", "--all",  dest = "all_compare", action="store_true", required= False, help="All files in supplied database will be compared")
compare_exclusive.add_argument("-s", "--samples",  dest = "samples_compare", metavar="sample_name[s]", nargs="+", required= False, help="Sample names supplied will be compared")

#extract_parser
extract_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database with the sample to remove")
extract_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the updated ddbb")

extract_parser.add_argument("-s", "--samples",  dest = "samples_extract", metavar="sample name[s]", nargs="+", required= True, help="Sample names supplied will be removed")


parser.add_argument("--version", action="version", version="%(prog)s 0.1")


arguments = parser.parse_args()





file = arguments.final_database
for i in arguments.samples_compare:
     print("samples to add are %s" %  i)

out_file = arguments.output_file

print("File is %s" % file)
print("Out File is %s" % out_file)