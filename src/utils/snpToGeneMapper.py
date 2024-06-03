import sys
from tqdm import tqdm

ENSAMBLE_PATH = "..\\data\\Ensembl\\PP_Ensembl_output.txt"
MAPPING_PATH = "..\\data\\snp_gene_mappings\\snpsToGene_PP_ensembl.txt"

if __name__ == "__main__":
    mapping = {}
    ensambleFile = open(ENSAMBLE_PATH, "r", encoding="utf-8")
    ensambleFile.readline()
    rows = ensambleFile.readlines()
    for row in tqdm(rows):
        line = row.strip().split("\t")
        snp = line[0]
        gene = line[5]
        if snp not in mapping:
            if gene == "-":
                mapping[snp] = "None"
            else:
                mapping[snp] = gene
        else:
            if gene == "-":
                continue
            else:
                if gene not in mapping[snp]:
                    mapping[snp] += "," + gene

    ensambleFile.close()

    print(mapping)
    
    with open(MAPPING_PATH, "w+", encoding="utf-8") as mappingFile:
        for snp in mapping:
            mappingFile.write(snp + " " + mapping[snp] + "\n")

    sys.exit()
            
