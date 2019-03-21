import pandas as pd
import os

def hamming_distance (pd_matrix):


    unique_values = pd.unique(pd_matrix[list(pd_matrix.keys())].values.ravel('K'))

    U = pd_matrix.eq(unique_values[0]).astype(int)
    H = U.dot(U.T)


    for unique_val in range(1,len(unique_values)):
        U = pd_matrix.eq(unique_values[unique_val]).astype(int)
        H = H.add(U.dot(U.T))

    return len(pd_matrix.columns) - H

#matrix_alleles = os.path.join('/media/bioinfo/NGS_Disk_Data','matrix_test.csv')
matrix_alleles = os.path.join('/home/smonzon/Downloads/','result_for_tree_diagram.tsv')

pd_matrix = pd.read_csv(matrix_alleles, sep='\t', header=0, index_col=0)

distance_matrix = hamming_distance (pd_matrix)
'''
unique_values = pd.unique(pd_matrix[list(pd_matrix.keys())].values.ravel('K'))
print ('hello')
U = pd_matrix.eq(unique_values[0]).astype(int)
H = U.dot(U.T)

#for uniq_val in unique_values[1:]:
#for unique_val in range (len(unique_values),1):

for unique_val in range(1,len(unique_values)):
    U = pd_matrix.eq(unique_values[unique_val]).astype(int)
    H = H.add(U.dot(U.T))
    print(unique_val)

'''
#print (distance_matrix)
distance_matrix.to_csv('distancia2.csv')
