import argparse


# parser = argparse.ArgumentParser()
# parser._action_groups.pop()
# required = parser.add_argument_group('required arguments')
# optional = parser.add_argument_group('optional arguments')
# required.add_argument('--required_arg')
# optional.add_argument('--optional_arg')
# return parser.parse_args()



# parser = argparse.ArgumentParser()
# optional = parser._action_groups.pop() # Edited this line
# required = parser.add_argument_group('required arguments')
# # remove this line: optional = parser...
# required.add_argument('--required_arg')
# optional.add_argument('--optional_arg')
# parser._action_groups.append(optional) # added this line
# return parser.parse_args(




parser = argparse.ArgumentParser(prog = 'VCF_DDTB.py', description= 'VCF_DDTB manages a custom snp database') 
#group = parser.add_mutually_exclusive_group()

#Positional arguments
parser.add_argument("-d", "--directory", help="Directory where vcf files are")
parser.add_argument("-o", "--out", help="Name of the new database", required=True)


parser.add_argument("-a", "--add", dest="new_file", help="File with variants to include in the database")
parser.add_argument("-c", "--compare", nargs='+' , type=str, metavar='SAMPLE NAMES', required= True, help="Compare sample names already in the database")
#action='append'

#default=42

parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument("-v","--verbose", action="store_true", help="increase output verbosity")
#choices=[0, 1, 2]

args = parser.parse_args()

if args.verbose:
    print("verbosity turned on")

file = args.new_file
for i in args.compare:
    print("samples to add are %s" %  i)

print("File is %s" % file)




https://github.com/BU-ISCIII/taranis/blob/master/taranis.py

TARANIS

#!/usr/bin/env python3
#import analyze_schema
import sys
import argparse
import allele_calling
import analyze_schema

def check_arg (args=None) :
    parser = argparse.ArgumentParser(prog = 'tara.py', formatter_class=argparse.RawDescriptionHelpFormatter, description= 'Taranis is a set of utilities related to cgMSLST ') 

    parser.add_argument('--version', action='version', version='%(prog)s 0.3.5')
    
    subparser = parser.add_subparsers(help = 'allele calling /interactive/schema are the available actions to execute taranis', dest = 'chosen_action')
    
    allele_calling_parser = subparser.add_parser('allele_calling', help = 'Allele calle way to downloads the  schema locus')
    allele_calling_parser.add_argument('-coregenedir', required= True, help = 'Directory where the core gene files are located ')
    allele_calling_parser.add_argument('-inputdir', required= True, help ='Directory where are located the sample fasta files')
    allele_calling_parser.add_argument('-outputdir', required= True, help = 'Directory where the result files will be stored')
    allele_calling_parser.add_argument('-cpus', required= False, help = 'Number of CPUS to be used in the program. Default is 1.', default = 1)
    allele_calling_parser.add_argument('-updateschema' , required=False, help = 'Create a new schema with the new locus found. Default is True.', default = True)
    allele_calling_parser.add_argument('-percentlength', required=False, help = 'Allowed length percentage to be considered as ASM or ALM. Outside of this limit it is considered as LNF Default is 20.', default = 20)
    
    evaluate_parser = subparser.add_parser('evaluate_schema', help = 'Evaluate the schema ')
    evaluate_parser.add_argument('-input_dir',required= True, help = 'Directory where are the schema files.')
    evaluate_parser.add_argument('-alt', required = False, help = 'Set to Yes if alternative start codon should be considered. Set to No to accept only ATG start codon', default = False)
    
    compare_parser = subparser.add_parser('compare_schema', help = 'Compare 2 schema')
    compare_parser.add_argument('-scheme1', help = 'Directory where are the schema files for the schema 1')
    compare_parser.add_argument('-scheme2', help = 'Directory where are the schema files for the schema 2')
    
    return parser.parse_args()



def processing_evaluate_schema (arguments) :
    print ('evaluate_schema')
    
    return True

def processing_compare_schema (arguments) :
    print ('compare_schema')
    return True

if __name__ == '__main__' :
    version = 'tara version 0.2.1'
    if len(sys.argv) == 1 :
        print( 'Mandatory parameters are missing to execute the program. \n ' ,'Usage: "tara.py  -help " for more information \n')
        exit (0)
    arguments = check_arg(sys.argv[1:])
    
    if arguments.chosen_action == 'allele_calling' :
        allele_calling.processing_allele_calling(arguments)
    elif arguments.chosen_action == 'evaluate_schema':
        analyze_schema.processing_evaluate_schema(arguments)
    elif arguments.chosen_action == 'compare_schema' :
        processing_compare_schema(arguments)
    print('completed') 