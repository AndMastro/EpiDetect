
# coding: utf-8

# # This module generates pickle file for the datasets

# ### import modules

import sys
import pickle
import csv
from tqdm import tqdm
from os import listdir
from os.path import isfile, join



DATASET_DIR_PATH = "../data/largeFiles/datasets/synthetic/epistatic_risk_model_MAF03_eta01_theta1_EDM-1/"
PICKLE_DIR_PATH = "../data/largeFiles/datasets/pickles/"



'''this function is employed to generate pickles with the format one-hot
    from datasetas created  using GAMETES. We keep the predictive SNPs as the lst two ones, since
    our MLP does not care about order. If using a different architecture, this should be revised.'''

def getDatasetFromGametes(saveMode=False, datasetPath=None, picklePath=None):
    
    dataset = {}
    numIndividuals = 0

    with open(datasetPath, mode='r') as data_file:
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
        with open(picklePath, 'wb+') as handle:
            pickle.dump(dataset, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Saved.")

if __name__ == "__main__":
    dataFiles = [f for f in listdir(DATASET_DIR_PATH) if isfile(join(DATASET_DIR_PATH, f))]
    for d in tqdm(dataFiles):
        fullPath = DATASET_DIR_PATH+d
        pPath = PICKLE_DIR_PATH + d[0:-3] + "p"
        getDatasetFromGametes(saveMode=True, datasetPath=fullPath, picklePath=pPath)
    print("All datsets generated.")
    sys.exit()