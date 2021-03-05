import csv
from os import sep
from tqdm import tqdm
import sys

no_van_snps = ["rs1041163", "rs1367117", "rs3741240", "rs1800590"]
FILE_NAME = "epistaticInteractions_NID_SBP_SUM_ordered.txt"
SAVE_FILE = "epistaticInteractions_NID_SBP_SUM_no_van_removed_ordered.txt"

saveFile = open(SAVE_FILE, mode="w+")

with open(FILE_NAME, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        for row in tqdm(csv_reader):
            snp1 = row[0]
            snp2 = row[1]
            if snp1 in no_van_snps or snp2 in no_van_snps:
                continue
            else:
                line = row[0] + "\t" + row[1] + "\t" + row[2] + "\t" + row[3] + "\t" + row[4] + "\n"
            saveFile.write(line) 