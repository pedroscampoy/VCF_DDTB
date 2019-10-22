#!/usr/bin/env python

import os
import pandas as pd
import argparse
import sys
import subprocess
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


#### OLD RECALIBRATING FUNCTIONS
########################################################################################################################################
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
##########################################################################################################################
#CREATED IN vcf_process.py

def import_VCF42_cohort_pandas(vcf_file, sep='\t'):
    """
    Script to read vcf 4.2 cohort/join called vcf handling header lines
    """
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            #print(next_line)
            next_line = f.readline()
    
    if first_line.endswith('VCFv4.2'):
        dataframe = pd.read_csv(vcf_file, sep=sep, skiprows=[header_lines], header=header_lines)
    else:
        print("This vcf file is not v4.2")
        sys.exit(1)

    return dataframe

def recheck_variant(format_sample):
    #GT:AD:DP:GQ:PGT:PID:PL:PS
    list_format = format_sample.split(":")
    gt = list_format[0]
    #gt0 = gt[0]
    #gt1 = gt[1]
    ad = list_format[1]
    ref = int(ad.split(',')[0])
    alt = max(int(x) for x in ad.split(',')[0:])
    
    if gt == "0/0":
        value = 0
    elif gt == "1/1":
        value = 1
    else:
        if gt == "./.":
            value = "!"
        elif "2" in gt:
            value = "!"
        elif (ref > alt):
            value = 0
        elif (alt > ref):
            value = 1
        else:
            value = "!"
            
    return value

def recheck_variant_mpileup(reference_file, position, sample, bam_folder):
    #Find reference name
    with open(reference_file) as f:
        reference = f.readline().split(" ")[0].strip(">").strip()
    #Identify correct bam
    for root, _, files in os.walk(bam_folder):
        for name in files:
            filename = os.path.join(root, name)
            if name.startswith(sample) and name.endswith(".bqsr.bam"):
                bam_file = filename
    #format position for mpileuo execution (NC_000962.3:632455-632455)
    position = reference + ":" + str(position) + "-" + str(position)
    
    #Execute command and retrieve output
    cmd = ["samtools", "mpileup", "-f", reference_file, "-aa", "-r", position, bam_file]
    text_mpileup = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, universal_newlines=True) 
    
    #Extract 5th column to find variants
    variant = text_mpileup.stdout.split()[4]
    var_list = list(variant)
    
    most_freq_var = max(set(var_list), key = var_list.count).upper()
        
    if most_freq_var == "." or most_freq_var == "," or most_freq_var == "*":
        return 0
    else:
        return 1

def identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, list_presence, bam_folder):
    """
    Replace nongenotyped ("!") with the most abundant genotype
    """
    #mode = max(set(list_presence), key = list_presence.count)
    
    count_ng = list_presence.count("!")
    sample_number = len(list_presence)
    
    if "!" not in list_presence:
        return list_presence
    elif count_ng/sample_number > 0.2:
        return 'delete'
    else:
        indices_ng = [i for i, x in enumerate(list_presence) if x == "!"]
        for index in indices_ng:
            #print(reference_file, row_position, sample_list_matrix[index], bam_folder)
            list_presence[index] = recheck_variant_mpileup(reference_file, row_position, sample_list_matrix[index], bam_folder)
        #new_list_presence = [mode if x == "!" else x for x in list_presence]
        return list_presence

def extract_recalibrate_params(pipeline_folder):
    for root, dirs, _ in os.walk(pipeline_folder):
        if root == pipeline_folder:
            for directory in dirs:
                subfolder = os.path.join(root, directory)
                if subfolder.endswith("/VCF"):
                    for file in os.listdir(subfolder):
                        if file.endswith("cohort.combined.hf.vcf"):
                            cohort_file = os.path.join(subfolder, file)
                            
                            with open(cohort_file, 'r') as f:
                                for line in f:
                                    if line.startswith("#"):
                                        if "--reference " in line:
                                            reference_file = line.split("--reference ")[1].strip().split(" ")[0].strip()
                                        
                            
                elif subfolder.endswith("/Bam"):
                    bam_folder = subfolder
                    
    return (cohort_file, bam_folder, reference_file)

def recalibrate_ddbb_vcf(snp_matrix_ddbb, vcf_cohort, bam_folder, reference_file):
    
    vcf_cohort = os.path.abspath(vcf_cohort)
    #snp_matrix_ddbb = os.path.abspath(snp_matrix_ddbb)
    
    df_matrix = snp_matrix_ddbb
    df_cohort = import_VCF42_cohort_pandas(vcf_cohort)
    
    sample_list_matrix = df_matrix.columns[3:]
    n_samples = len(sample_list_matrix)
    #sample_list_cohort = df_cohort.columns.tolist()[9:]
    
    
    list_index_dropped = []
    #Iterate over non unanimous positions 
    for index, data_row in df_matrix[df_matrix.N < n_samples].iloc[:,3:].iterrows():
        #Extract its position
        row_position = int(df_matrix.loc[index,"Position"])
        #print(data_row.values)
        #Use enumerate to retrieve column index (column ondex + 3)
        presence_row = [recheck_variant(df_cohort.loc[df_cohort.POS == row_position, df_matrix.columns[n + 3]].item()) \
                           for n,x in enumerate(data_row)]
        #print(presence_row, row_position)
        #Resolve non genotyped using gvcf files
        new_presence_row = identify_nongenotyped_mpileup(reference_file, row_position, sample_list_matrix, presence_row, bam_folder)
        
        #find positions with 20% of nongenotyped and delete them OR
        #reasign positions without nongenotyped positions 
        if new_presence_row == 'delete':
            list_index_dropped.append(index)
        else:
            df_matrix.iloc[index, 3:] = new_presence_row
            df_matrix.loc[index, 'N'] = sum(new_presence_row)
        #print(new_presence_row)
        #print("\n")
    #Remove all rows at once to avoid interfering with index during for loop
    df_matrix.drop(index=list_index_dropped, axis=0, inplace=True)
    
    return df_matrix


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
        args.recalibrate = os.path.abspath(args.recalibrate)
        if os.path.exists(args.recalibrate):
            recalibrate_params = extract_recalibrate_params(args.recalibrate)
            print("\n" + MAGENTA + "Recalibration selected" + END_FORMATTING)
            output_file = (".").join(output_file.split(".")[:-1]) + ".revised.tsv"
            final_ddbb_revised = recalibrate_ddbb_vcf(final_ddbb, recalibrate_params[0], recalibrate_params[1], recalibrate_params[2])
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