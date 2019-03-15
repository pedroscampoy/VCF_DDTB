import argparse



parser = argparse.ArgumentParser(prog = 'VCF_DDTB.py', formatter_class=argparse.RawDescriptionHelpFormatter, description= 'VCF_DDTB manages a custom snp database') 

#Positional arguments
parser.add_argument("-n", "--new", help="File with variants to include in the database")
parser.add_argument("-v","--verbose", action="store_true", help="increase output verbosity")
#choices=[0, 1, 2]

args = parser.parse_args()

if args.verbose:
    print("verbosity turned on")

file = args.new

print("File is %s" % file)