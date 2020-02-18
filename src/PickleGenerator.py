
# coding: utf-8

# # This module generates pickle file for the datasets

# ### import modules

import sys
import pickle
import csv
from tqdm import tqdm



DATASET_DIR_PATH = "../data/largeFiles/datasets/synthetic/"
DATASET_NAME = "prova_EDM-1/prova_EDM-1_1.txt"
DATASET_PATH = DATASET_DIR_PATH + DATASET_NAME
PICKLE_PATH = "../data/largeFiles/datasets/pickles/prova.p"


'''this function is emplyed to generate pickles with the format one-hot
    from datasetas created  using GAMETES. We keep the predictive SNPs as the lst two ones, since
    our MLP does not care about order. If using a different architecture, this should be revised.'''

def getDatasetFromGametes(saveMode=False):
    
    dataset = {}
    numIndividuals = 0

    with open(DATASET_PATH, mode='r') as data_file:
        csv_reader = csv.reader(data_file, delimiter = "\t")
        next(csv_reader) #skip first line since it's header

        print("Generating dataset...")

        for row in tqdm(csv_reader):
            phenotype = int(row[-1]) #class label is the last elem of the arrayù
            genotype = [] #genotype must be encoder as one-hot encoding
            for snp in row[0:-1]:
                if snp == "0":
                    genotype += [1,0,0]
                elif snp == "1":
                    genotype += [0,1,0]
                elif snp == "2":
                    genotype += [0,0,1]
                else:
                    genotype += [0,0,0] #just for robustness, there are no missing data in synthatic datasets.
                
            dataset[numIndividuals] = [genotype, phenotype]
            numIndividuals += 1

    if saveMode == True:
        print("Saving pickle...")
        with open(PICKLE_PATH, 'wb+') as handle:
            pickle.dump(dataset, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Saved.")

if __name__ == "__main__":
    getDatasetFromGametes(saveMode=True)
    sys.exit()