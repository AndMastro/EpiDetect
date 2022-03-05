
# coding: utf-8

# # This module generates pickle file for the datasets

# ### import modules

import sys
import pickle
import csv
from tqdm import tqdm
from os import listdir
from os.path import isfile, join


DATASET_TYPE = "marginal_effect"
DATASET_DIR_PATH = "../data/largeFiles/datasets/GAMETES/" + DATASET_TYPE + "/epistatic_plus_ME_risk_model_MAF01_eta01_theta1_lambda03_1_2_EDM-1/"
TEXT_DIR_PATH = "../data/largeFiles/datasets/text/GAMETES/" + DATASET_TYPE + "/"



'''this function is employed to generate text with the format for BOOST
    from datasetas created  using GAMETES. We keep the predictive SNPs as the last two ones.
    Phenotype is the first elem in row.'''

if __name__ == "__main__":
    dataFiles = [f for f in listdir(DATASET_DIR_PATH) if isfile(join(DATASET_DIR_PATH, f))]
    for d in tqdm(dataFiles):
        fullPath = DATASET_DIR_PATH+d
        
        savePath = TEXT_DIR_PATH + d[0:-4] + "_BOOST.txt"
        saveFile = open(savePath, "w+")
        firstLine = False
        with open(fullPath, "r") as dataFile:
            print("\nReading file " + d)
            for line in tqdm(dataFile):
                
                if not firstLine:
                    firstLine = True
                    continue
            
                row = line.strip().split("\t")
                phenotype = row[-1]
                genotype = row[:-1]
                numSNPs = len(genotype)
                added = 0
                saveFile.write(phenotype + " ") #phenotype is the first elem in row for BOOST (oppositely eith MDR, where it is the last)
                for snp in genotype:
                    added += 1
                    saveFile.write(snp)
                    if added == numSNPs:
                        saveFile.write("\n")
                    else:
                        saveFile.write(" ")
        

    print("All datasets generated.")
    sys.exit()