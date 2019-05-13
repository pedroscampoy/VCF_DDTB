#!/home/laura/env36/bin/python

import os
import pandas as pd
import argparse
import sys
from misc import check_file_exists, import_to_pandas, extract_sample_snp_final


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


def ddtb_new(args):
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
        if not filename.startswith('.') and filename.endswith(".snp.final"):
            print("\nThe file is: %s" % filename)
            
            all_samples = all_samples + 1
            positions_shared = []
            positions_added = []
            
            sample = extract_sample_snp_final(filename) #Manage sample name
            
            file = os.path.join(directory, filename) #Whole file path
            check_file_exists(file) #Manage file[s]. Check if file exist and is greater than 0

            new_sample = import_to_pandas(file) #Import files in annotated vcf format

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
                for position in new_sample.iloc[:,0].unique(): #extract first column in file
                    
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
                
            #Create small report with basic count
            #####################################
            
    print("\n" + GREEN + "Position check Finished" + END_FORMATTING)
    print(GREEN + "Added " + str(new_samples) + " samples out of " + str(all_samples) + END_FORMATTING + "\n")
    #Sort positions and save new DDBB to TSV
    final_ddbb['N'] = final_ddbb['N'].astype(int)
    #pd.set_option('display.precision', 0)
    #pd.reset_option('^display.', silent=True) #Reset options in case I mess up

    final_ddbb = final_ddbb.fillna(0).sort_values("Position") #final_ddbb = final_ddbb["Position"].astype(int)
    print("Final database now contains %s rows and %s columns" % final_ddbb.shape)

    final_ddbb.to_csv(output_file, sep='\t', index=False)