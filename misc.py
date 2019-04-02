import pandas as pd
import os
import sys

def check_file_exists(file_name):
    """
    Check file exist, if not program exit.
    """
    if not os.path.isfile(file_name):
        sys.exit(1)

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