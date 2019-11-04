#!/usr/bin/env python

import os
import pandas as pd
from sklearn.metrics import jaccard_score, pairwise_distances
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd #pdist
import PyQt5
from PyQt5 import QtGui
import ete3
from ete3 import Tree, TreeStyle

from misc import check_file_exists, import_to_pandas

"""
TODO:   confirm 'average' method as the best option
        confirm euclidean 'metric' as the most appropiate for SNP distance (hamming used so far)
"""

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


def compare_snp_columns(sample1, sample2, df):
    jaccard_similarity = jaccard_score(df[sample1], df[sample2], average='binary') #similarities between colums
    hamming_similarity = 1 - jaccard_similarity #disagreements between colums
    snp_distance = int(hamming_similarity * (len(df.index)+1))
    return snp_distance

def snp_distance_pairwise(dataframe, output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    with open(output_file, "a") as f:
        for sample1 in dataframe.iloc[:,3:].columns: #remove first 3 colums
            for sample2 in dataframe.iloc[:,3:].columns:
                if sample1 != sample2:
                    snp_distance = compare_snp_columns(sample1, sample2, dataframe)
                    line_distance = "%s\t%s\t%s\n" % (sample1, sample2, snp_distance)
                    f.write(line_distance)

def snp_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index), index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    snp_distance_df = snp_distance_df.astype(int)
    snp_distance_df.to_csv(output_file, sep='\t', index=True)

def hamming_distance_matrix(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(dataframe_only_samples.T, metric = "hamming") #dataframe.T means transposed
    hamming_distance_df = pd.DataFrame(hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns) #Add index
    hamming_distance_df.to_csv(output_file, sep='\t', index=True)

def clustermap_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    sns.clustermap(dataframe_only_samples, annot=False, cmap="YlGnBu", figsize=(13, 13))
    plt.savefig(output_file, format="png")

def dendogram_dataframe(dataframe, output_file):
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average') #method='single'

    plt.rcParams['lines.linewidth'] = 8 #Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10 #Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30}) #Increase x tick label size
    #plt.tick_params(labelsize=30)
    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending', show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    
    plt.savefig(output_file, format="png")

# Convert dendrogram to Newick
def linkage_to_newick(dataframe, output_file):
    """
    Thanks to https://github.com/biocore/scikit-bio/issues/1579
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """
    dataframe_only_samples = dataframe.set_index(dataframe['Position'].astype(int)).drop(['Position','N','Samples'], axis=1) #extract three first colums and use 'Position' as index
    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average')

    tree = shc.to_tree(Z, False)
    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            #print("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            #print(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)



def ddtb_compare(args):

    database_file = os.path.abspath(args.final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    if args.output_file:
        output_file = os.path.abspath(args.output_file)
        output_path = output_file.split(".")[0]
    else:
        output_path = database_file.split(".")[0]

    print("Output path is: " + output_path)

    if args.all_compare:
        print(BLUE + BOLD + "Comparing all samples in " + database_file + END_FORMATTING)
        prior_pairwise = datetime.datetime.now()

        #Calculate pairwise snp distance for all and save file
        print(CYAN + "Pairwise distance" + END_FORMATTING)
        pairwise_file = output_path + ".snp.pairwise.tsv"
        snp_distance_pairwise(presence_ddbb, pairwise_file)
        after_pairwise = datetime.datetime.now()
        print("Done with pairwise in: %s" % (after_pairwise - prior_pairwise))

        #Calculate snp distance for all and save file
        print(CYAN + "SNP distance" + END_FORMATTING)
        snp_dist_file = output_path + ".snp.tsv"
        snp_distance_matrix(presence_ddbb, snp_dist_file)

        #Calculate hamming distance for all and save file
        print(CYAN + "Hamming distance" + END_FORMATTING)
        hmm_dist_file = output_path + ".hamming.tsv"
        hamming_distance_matrix(presence_ddbb, hmm_dist_file)
        """
        #Represent pairwise snp distance for all and save file
        print(CYAN + "Drawing distance" + END_FORMATTING)
        prior_represent = datetime.datetime.now()
        png_dist_file = output_path + ".snp.distance.png"
        #clustermap_dataframe(presence_ddbb, png_dist_file)
        after_represent = datetime.datetime.now()
        print("Done with distance drawing in: %s" % (after_represent - prior_represent))
        """
        #Represent dendrogram snp distance for all and save file
        print(CYAN + "Drawing dendrogram" + END_FORMATTING)
        png_dend_file = output_path + ".snp.dendrogram.png"
        dendogram_dataframe(presence_ddbb, png_dend_file)


        #Output a Newick file distance for all and save file
        print(CYAN + "Newick dendrogram" + END_FORMATTING)
        newick_file = output_path + ".nwk"
        linkage_to_newick(presence_ddbb, newick_file)


    else:
        print("sample mode is not implemented")







    """
    #compare_parser
     compare_parser.add_argument("-d", "--database",  dest = "final_database", required= True, metavar="TB_database", help="REQUIRED: csv file with the database to be enriched/consulted")
     compare_parser.add_argument("-o", "--outputfile",  dest = "output_file", required= True, metavar="filename", help="REQUIRED: file name, including PATH, of the matrix comparison")

     compare_exclusive = compare_parser.add_mutually_exclusive_group()

     compare_exclusive.add_argument("-a", "--all",  dest = "all_compare", action="store_true", required= False, help="All files in supplied database will be compared")
     compare_exclusive.add_argument("-s", "--samples",  dest = "samples_compare", metavar="sample_name[s]", nargs="+", required= False, help="Sample names supplied will be compared")

    """