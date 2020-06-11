from snpReader import geneRead
import csv
from tqdm import tqdm

if __name__ == "__main__":
    with open("../data/largeFiles/PP_2_measures.csv", mode='r') as csv_file:
         csv_reader = csv.reader(csv_file)
         for row in tqdm(csv_reader):
            if (row[0] == "-1"):
                print(row)
                break
