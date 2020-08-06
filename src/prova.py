from snpReader import geneRead
import csv
from tqdm import tqdm
import pickle
if __name__ == "__main__":

    data = None
    with open("../data/largeFiles/pickles/BP/SNPS_PP_AVG.p", "rb") as dataFile:
        data = pickle.load(dataFile)

    INDS = list(data.keys())
    fv = data[INDS[0]][0]
    print(len(fv)/3)

    '''with open("../data/largeFiles/pickles/BP", mode='r') as csv_file:
         csv_reader = csv.reader(csv_file)
         for row in tqdm(csv_reader):
            if (row[0] == "-1"):
                print(row)
                break'''
