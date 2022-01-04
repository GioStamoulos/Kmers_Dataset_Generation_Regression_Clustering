from itertools import product
import pandas as pd
from collections import defaultdict
import os
import csv

class kmers_dataset:

    def __init__(self, k_mers, output_folder, genome_folder):
        self.k_mers = k_mers
        self.output_folder = output_folder
        self.genome_folder = genome_folder
        
    def create_csv_headers(self):
        letters =['A', 'C', 'G', 'T']
        header_list = [''.join(i) for i in product(letters, repeat = self.k_mers)]
        header_list.insert(0, 'Genome_Name')
        return(header_list)


    def create_dataset(self, header_list):
        self.header_list = header_list
        # initializing dictationary
        d_con = defaultdict(list)
        for entry in os.scandir(self.output_folder):
            genome_name =  str(os.path.basename(entry))[:-4]
            d = pd.read_csv(entry, header= None, sep='\t', index_col=0).to_dict() 
            d_temp = d[1]
            d_temp['Genome_Name'] = genome_name
            for i in self.header_list:
                try:
                    d_con[i].append(d_temp[i])
                except:
                    d_con[i].append(0)
        #create dataframe from dictationary
        df = pd.DataFrame.from_dict(d_con)
        #create kmers dataset - csv file 
        df.to_csv('{0}/{1}/kmers_dataset.csv'.format(self.genome_folder, self.k_mers), index=False)
        print(d_con)                         