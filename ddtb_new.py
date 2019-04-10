#!/home/laura/env36/bin/python

import os
import pandas as pd
import argparse
from misc import check_file_exists, import_to_pandas

def blank_database():
    new_pandas_ddtb = pd.DataFrame(columns=['Position','N', 'Samples'])
    return new_pandas_ddtb


def ddtb_new(args):
    directory = args.folder
    #cwd = os.getcwd()

    final_ddbb = blank_database()

    print("Previous final database contains %s rows and %s columns\n" % final_ddbb.shape)

    print("The directory selected is: %s" % directory)


    for filename in os.listdir(directory):
        if not filename.startswith('.'):
            print("The file is: %s" % filename)
            
            positions_shared = []
            positions_added = []
            
            #Manage sample name. Split by "_" and take the first word on the left
            ###################
            sample = filename.split("_")[0]        
            
            #Manage file[s]. Check if file exist and is greater than 0
            ###############
            file = os.path.join(directory, filename) #Whole file path
            file_info = os.stat(file) #Retrieve the file info to check if has size > 0


            if os.path.isfile(file) and file_info.st_size > 0:
                
                pass

            else:
                print("ERROR: Your file %s does not exist or is empty" % filename)

            #Import files in annotated vcf format
            #####################################

            #new_sample = pd.read_csv(file, sep='\t', skiprows=[0], header=None)
            new_sample = import_to_pandas(file)

            #Handle each new_sample
            #######################
            
            print("This file contains %s SNPs" % len(new_sample.index))
            
            #Check if sample exist
            ######################
            
            if sample not in final_ddbb.columns.tolist():
                print("The sample %s is NOT in final ddbb. Adding sample" % sample)
                #print("Adding new sample %s to final ddbb" % sample)
                
                #extrac the number of columns to insert a new one
                new_colum_index = len(final_ddbb.columns)
                #final_ddbb[sample] = sample #adds a new column but fills all blanks with the value sample
                
                #add a new column with defauls values = 0
                final_ddbb.insert(new_colum_index, sample, 0)
                
                #Check if position exist
                ########################
                
                print("Checking positions(SNPs) in sample %s" % sample)
                for position in new_sample.iloc[:,0].unique(): #extract first column in file

                    #print(type(position))
                    
                    if position not in final_ddbb["Position"].values:
                    
                        positions_added.append(position)
                        
                        new_row = len(final_ddbb.index)
                        final_ddbb.loc[new_row,'Position'] = position
                        final_ddbb.loc[new_row,'Samples'] = sample
                        final_ddbb.loc[new_row,'N'] = int(1)
                        final_ddbb.loc[new_row,sample] = str(1)
                        
                    else:
                        
                        positions_shared.append(position)
                        
                        #Check whether the column matches the value and retrieve the first position [0]
                        #of the object index generated
                        index_position = final_ddbb.index[final_ddbb["Position"] == position][0]
                        
                                                    
                        number_samples_with_position = final_ddbb.loc[index_position,'N']
                        names_samples_with_position = final_ddbb.loc[index_position,'Samples']
                        new_names_samples = names_samples_with_position + "," + sample
                        
                        #print("NEW NAME = %s" % new_names_samples)
                        
                        final_ddbb.loc[index_position,'N'] = number_samples_with_position + 1 #add 1 to the number of samples
                        final_ddbb.loc[index_position,'Samples'] = new_names_samples
                        final_ddbb.loc[index_position,sample] = int(1) #Add "1" in cell with correct position vs sample
            
            else:
                
                print("This sample ALREADY exist, have a look")
                    
            #Create small report with basic count
            #####################################
            
            print("\nSAMPLE:\t%s\nTOTAL Variants:\t%s\n\Shared Variants:\t%s\nNew Variants:\t%s\n" % (sample, len(new_sample.index), len(positions_shared), len(positions_added)))

    #final_ddbb = final_ddbb["Position"].astype(int)
    pd.set_option('display.precision', 0)
    #pd.reset_option('^display.', silent=True) #Reset options in case I mess up

    final_ddbb = final_ddbb.fillna(0).sort_values("Position")
        
    print("Final database now contains %s rows and %s columns" % final_ddbb.shape)
    #print(final_ddbb)

    final_ddbb.to_csv(args.output_file, sep='\t', index=False)