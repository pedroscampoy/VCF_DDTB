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
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=0)
    
    return dataframe
