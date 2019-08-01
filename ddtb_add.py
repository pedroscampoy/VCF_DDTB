#!/usr/bin/env python

import os
import pandas as pd
import argparse
import sys
from misc import check_file_exists, import_to_pandas, extract_sample_snp_final, import_VCF4_to_pandas

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


def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb


def retrieve_tabs(sample_list, folder_tab):
    dict_tab_files = {}
    for sample in sample_list:
        for root, _, files in os.walk(folder_tab):
            for name in files:
                filename = os.path.join(root, name)
                if (sample in filename) and filename.endswith(".tab"):
                    dict_tab_files[sample] = pd.read_csv(filename, sep="\t",header=0)
    if len(sample_list) == len(dict_tab_files.keys()):
        return dict_tab_files
    else:
        print('Some tab files are missing, please make sure all samples are represented')
        sys.exit(1)

def variant_is_present(vcf_df,position):
    if position in vcf_df.POS.tolist():
        index_position = vcf_df.index[vcf_df.POS == position][0]
        if vcf_df.loc[index_position,'ALT_AD'] > vcf_df.loc[index_position,'REF_AD']:
            return "1"
        else:
            return "0"
    else:
        return "0"

def recallibrate_ddbb(snp_matrix_ddbb, folder_tab):
    
    folder_tab = os.path.abspath(folder_tab)
    df = snp_matrix_ddbb
    #df = pd.read_csv(snp_matrix_ddbb, sep="\t",header=0)
    sample_list = df.columns[3:]
    n_samples = len(sample_list)
    
    all_samples_tab = retrieve_tabs(sample_list,folder_tab)
    
    for sample in sample_list:
        for index, _ in df[df.N < n_samples].iterrows():
            previous_snp = df.loc[index,sample]
            #df.loc[index,sample] coordinates to replace
            #all_samples_tab[sample] Dictionary with tab dataframes with all samples and variants
            #df.iloc[index,0] First columnn = POSITION to check with function and return new 0 or 1
            post_snp = int(variant_is_present(all_samples_tab[sample], df.iloc[index,0]))
            df.loc[index,sample] = post_snp
            #df.loc[index,sample] = variant_is_present(all_samples_tab[sample], df.iloc[index,0])
            #Substitute previous count (N) and list of samples
            if previous_snp == 0 and post_snp == 1:
                #Reassign number of samples (colimn with index 1)
                df.iloc[index,1] = df.iloc[index,1] + 1
                #Reassign list of samples (colimn with index 2)
                df.iloc[index,2] = df.iloc[index,2] + "," + sample
    return df

def ddtb_add(args):
    directory = os.path.abspath(args.folder)
    output_file = os.path.abspath(args.output_file)

    #Select NEW vs UPDATE
    if args.subtask == 'new' :
        final_ddbb = blank_database()
    elif args.subtask == 'update':
        update_database = os.path.abspath(args.update_database)
        if update_database == output_file:
            print(RED + "ERROR: " + END_FORMATTING + BOLD + "Pick a diferent name for the output database" + END_FORMATTING)
            sys.exit(1)
        else:
            final_ddbb = import_to_pandas(update_database, header=True)
    #Make sure output exist to force change name
    if os.path.isfile(output_file):
        print(YELLOW + "ERROR: " + BOLD + "output database EXIST, choose a different name or manually delete" + END_FORMATTING)
        sys.exit(1)

    
    print("Previous final database contains %s rows and %s columns\n" % final_ddbb.shape)
    print("The directory selected is: %s" % directory)
    

    all_samples = 0
    new_samples = 0
    for filename in os.listdir(directory):
        if not filename.startswith('.') and filename.endswith(args.suffix):
            print("\nThe file is: %s" % filename)
            
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            
            sample = filename.split(".")[0] #Manage sample name
            
            file = os.path.join(directory, filename) #Whole file path
            check_file_exists(file) #Manage file[s]. Check if file exist and is greater than 0

            new_sample = import_VCF4_to_pandas(file) #Import files in annotated vcf format

            #Handle each new_sample
            #print("This file contains %s SNPs" % len(new_sample.index))
            
            #Check if sample exist
            ######################
            if sample not in final_ddbb.columns.tolist():
                print("Adding new sample %s to %s" % (sample, os.path.basename(args.output_file)))
                new_samples = new_samples + 1
                new_colum_index = len(final_ddbb.columns) #extract the number of columns to insert a new one
                #final_ddbb[sample] = sample #adds a new column but fills all blanks with the value sample
                final_ddbb.insert(new_colum_index, sample, 0) #add a new column with defauls values = 0
                
                #Check if position exist
                ########################
                for position in new_sample['POS'].unique(): #extract first column in file
                    
                    if position not in final_ddbb["Position"].values:
                        positions_added.append(position) #Count new positions for stats
                        
                        new_row = len(final_ddbb.index)
                        final_ddbb.loc[new_row,'Position'] = position
                        final_ddbb.loc[new_row,'Samples'] = sample
                        final_ddbb.loc[new_row,'N'] = int(1)
                        final_ddbb.loc[new_row,sample] = str(1)
                        
                    else:
                        positions_shared.append(position) #Count shared positions for stats
                        
                        #Check whether the column matches the value and retrieve the first position [0]
                        #of the object index generated
                        index_position = final_ddbb.index[final_ddbb["Position"] == position][0]
                        #Add sample to corresponding cell [position, samples]
                        number_samples_with_position = final_ddbb.loc[index_position,'N']
                        names_samples_with_position = final_ddbb.loc[index_position,'Samples']
                        new_names_samples = names_samples_with_position + "," + sample
                        #Sum 1 to the numbes of samples containing the position
                        final_ddbb.loc[index_position,'N'] = number_samples_with_position + 1
                        final_ddbb.loc[index_position,'Samples'] = new_names_samples
                        final_ddbb.loc[index_position,sample] = str(1) #Add "1" in cell with correct position vs sample (indicate present)

                print("\nSAMPLE:\t%s\nTOTAL Variants:\t%s\nShared Variants:\t%s\nNew Variants:\t%s\n"
                % (sample, len(new_sample.index), len(positions_shared), len(positions_added)))
            else:
                print(YELLOW + "The sample " + sample + " ALREADY exist" + END_FORMATTING)

    final_ddbb = final_ddbb.fillna(0).sort_values("Position") #final_ddbb = final_ddbb["Position"].astype(int)
    final_ddbb['N'] = final_ddbb['N'].astype(int)
    final_ddbb = final_ddbb.reset_index(drop=True)

    print("Final database now contains %s rows and %s columns" % final_ddbb.shape)
    if args.recalibrate == False:
        final_ddbb.to_csv(output_file, sep='\t', index=False)
    else:
        if os.path.exists(args.recalibrate):
            print("\n" + MAGENTA + "Recalibration selected" + END_FORMATTING)
            output_file = (".").join(output_file.split(".")[:-1]) + ".revised.tsv"
            final_ddbb_revised = recallibrate_ddbb(final_ddbb, args.recalibrate)
            final_ddbb_revised.to_csv(output_file, sep='\t', index=False)
        else:
            print("The directory supplied for recalculation does not exixt")
            sys.exit(1)

    #Create small report with basic count
    #####################################
            
    print("\n" + GREEN + "Position check Finished" + END_FORMATTING)
    print(GREEN + "Added " + str(new_samples) + " samples out of " + str(all_samples) + END_FORMATTING + "\n")
    
    #pd.set_option('display.precision', 0)
    #pd.reset_option('^display.', silent=True) #Reset options in case I mess up