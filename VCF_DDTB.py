import os
import pandas as pd


directory = "test/test"
print(directory)
for filename in os.listdir(directory):

    print(filename)

    file = os.path.join(directory, filename)
    file_info = os.stat(file)

    print(file)

    if os.path.isfile(file) and file_info.st_size > 0:
        #continue
        print(filename)
    else:
        print("WARNING: Your file %s does not exist or is empty" % filename)

    #import files

    ddbb_tb = pd.read_csv(file, sep='\t', header=0)
    
    #new_sample
    
    print("This file contains %s rows and %s columns" % ddbb_tb.shape)

    #print(ddbb_tb.columns)

    #for i in ddbb_tb.iloc[:,0]:
    #for i in ddbb_tb["Position"]:
     #   print("Hola %s" % i)

    #check new variables

    #Add samples and increase N



