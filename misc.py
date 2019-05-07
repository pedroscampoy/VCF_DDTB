import pandas as pd
import os
import sys
import re

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def check_file_exists(file_name):
    """
    Check file exist and is not 0 Kb, if not program exit.
    """
    file_info = os.stat(file_name) #Retrieve the file info to check if has size > 0

    if not os.path.isfile(file_name) or file_info.st_size > 0:
        print(RED + BOLD + "File: %s not found or empty\n" % file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def import_to_pandas(file_table, header=False, sep='\t'):
    if header == False:
        #exclude first line, exclusive for vcf outputted by PipelineTB
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        #Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)
    
    return dataframe

def remove_if_exixt(file):
    if os.path.exists(file):
        os.remove(file)

def extract_sample_snp_final(snp_final):
    """
    Extract sample from snp.final file.
    """

    sample_name_R = snp_final.split("_")[0]
  
    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    
    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.rstrip(match)
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.rstrip(match)
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip(match)
    else:
        sample_name = sample_name_R

    return sample_name