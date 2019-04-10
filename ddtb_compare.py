#!/home/laura/env36/bin/python

import os
import pandas as pd
from sklearn.metrics import jaccard_similarity_score, pairwise_distances
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from misc import check_file_exists, import_to_pandas

def compare_jaccard_columns(sample1, sample2, df):
    jaccard_similarity = jaccard_similarity_score(df[sample1], df[sample2])
    return jaccard_similarity

def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = jaccard_similarity_score(df[sample1], df[sample2]) #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns: #remove first 3 colums
            for sample2 in dataframe.iloc[:,3:].columns:
                snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                f.write(line_distance) 
                #print(sample1, sample2, snp_distance)

def snp_distance_matrix(dataframe):
    data_hamming = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(data_hamming.T, metric = "hamming") #dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(data_hamming.index), index=data_hamming.columns, columns=data_hamming.columns) #Add index
    snp_distance_df = snp_distance.astype(int)
    return snp_distance_df



def ddtb_compare(args):
    if args.all_compare:
        snp_distance_matrix(args.databse)
    print("You choose COMPARE - This option is not yet implemented - sorry")







    """
    #compare_parser
     compare_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
     compare_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the matrix comparison")

     compare_exclusive = compare_parser.add_mutually_exclusive_group()

     compare_exclusive.add_argument("-a", "--all",  dest = "all_compare", action="store_true", required= False, help="All files in supplied database will be compared")
     compare_exclusive.add_argument("-s", "--samples",  dest = "samples_compare", metavar="sample_name[s]", nargs="+", required= False, help="Sample names supplied will be compared")

    """