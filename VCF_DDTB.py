import os
#from pathlib import Path
import pandas as pd


directory = "test"

for filename in os.listdir(directory):

    file = os.path.join(directory, filename)
    file_info = os.stat(file)

    if os.path.isfile(file) and file_info.st_size > 0:
        continue
        # print(file)
    else:
        print("Your file %s does not exist or is empty" % filename)


# else:
    #    print("NO")
    # if filename.endswith:
        # print(os.path.join(directory, filename))
    #    continue
    # else:
    # continue
