import os
import pandas as pd
import argparse
from misc import check_file_exists, import_to_pandas

file_for_comparing = 'test_100.csv'

df_for_comparing = import_to_pandas(file_for_comparing, header=True)

print(df_for_comparing.head())