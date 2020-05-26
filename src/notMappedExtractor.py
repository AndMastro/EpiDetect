import sys
from tqdm import tqdm

MAPPING_PATH = "..\\data\\snp_gene_mappings\\snpsToGene_DBP_ensembl.txt"
MISSING_PATH = "..\\data\\snp_gene_mappings\\notMapped_DBP_ensembl.txt"

if __name__ == "__main__":
    notMapped = []
    mapFile = open(MAPPING_PATH, "r", encoding="utf-8")
    rows = mapFile.readlines()
    for row in tqdm(rows):
        line = row.strip().split(" ")
        snp = line[0]
        gene = line[1]
        if gene == "None":
            notMapped.append(snp)

    mapFile.close()

    print(notMapped)
    
    with open(MISSING_PATH, "w+", encoding="utf-8") as missingFile:
        for snp in notMapped:
            missingFile.write(snp + "\n")

    sys.exit()
            
