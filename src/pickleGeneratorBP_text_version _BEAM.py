import sys
import csv
import pandas
import pickle
from tqdm import tqdm
from snpReader import snpRead

#SBP data
IND_GENO_PATH_SQ = "..\\data\\largeFiles\\allChrom_SBP_recoded12.csv"
IND_SBP = "..\\data\\largeFiles\\SBP_2_measures.csv" 
TEXT_PATH = "..\\data\\largeFiles\\datasets\\text\\BP\\SNPS_SBP_AVG_BEAM.txt"
SNP_FILE = "..\\data\\allChrom_SBP.bim"

#DBP data
'''IND_GENO_PATH_SQ = "..\\data\\largeFiles\\allChrom_DBP_recoded12.csv"
IND_SBP = "..\\data\\largeFiles\\DBP_2_measures.csv" 
TEXT_PATH = "..\\data\\largeFiles\\datasets\\text\\BP\\SNPS_DBP_AVG_BEAM.txt"
SNP_FILE = "..\\data\\allChrom_DBP.bim"'''

#PP data
'''IND_GENO_PATH_SQ = "..\\data\\largeFiles\\allChrom_PP_recoded12.csv"
IND_SBP = "..\\data\\largeFiles\\PP_2_measures.csv" 
TEXT_PATH = "..\\data\\largeFiles\\datasets\\text\\BP\\SNPS_PP_AVG_BEAM.txt"
SNP_FILE = "..\\data\\allChrom_PP.bim"'''


if __name__ == "__main__":

    
    genoType = {}
    
    print("Loading individuals genotypes...")
    with open(IND_GENO_PATH_SQ, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in tqdm(csv_reader):
            indID = row[0]
            genoVector = []
            dim = len(row)
            for i in range(1,dim-1, 2):
                if row[i] == "1" and row[i+1] == "2":
                    genoVector += [1]
                elif row[i] == row[i+1] == "1":
                    genoVector += [0]
                elif row[i] == row[i+1] == "2":
                    genoVector += [2]
                else:
                    genoVector += [-1] #BEAM wants negative numbers for missing alleles
            genoType[indID] = genoVector

    print(genoType["1000018"])

    completeDataset= {}
    print("Generating full dataset with genotype, phenotypes and BP values...")
    with open(IND_SBP, mode='r') as csv_file:
         csv_reader = csv.reader(csv_file)
         for row in tqdm(csv_reader):
            if row[0] in genoType:
                indID = row[0]
                genoVec = genoType[indID]
                phenoVec = []
                #phenoVec = row[1:-1] used only for phenotype
                sbpValue_0 = row[1]
                sbpValue_1 = row[2]
                sbpValue = 0
                if sbpValue_1 != "0":
                    sbpValue = (float(sbpValue_0) + float(sbpValue_1)) / 2.0
                else:
                    sbpValue = sbpValue_0
               
                #we look in the phenoVec list for compatibility, eve if we do not use pheno data
                completeDataset[indID] = [genoVec, [float(elem) for elem in phenoVec], float(sbpValue)]
    
    print(completeDataset["1596781"])

    for ind in completeDataset:
        print(completeDataset[ind][0])
        break
    
    print("Reading snps...")
    sbpSNPs = snpRead(SNP_FILE)

    
    print("Dataset dict generated. Saving as text...")

    with open(TEXT_PATH, 'w+') as saveFile:
        numSNPs = len(sbpSNPs)
        added = 0
        for snp in sbpSNPs:
            saveFile.write(sbpSNPs[snp])
            added += 1
            saveFile.write("\t")
            if added == numSNPs:
                saveFile.write("Outcome" + "\n")

        for ind in tqdm(completeDataset):
            geno = completeDataset[ind][0]
            BP = completeDataset[ind][2]
            for elem in geno:
                saveFile.write(str(elem) + "\t")
            saveFile.write(str(BP) + "\n")

    print("File saved.")

    sys.exit()

   