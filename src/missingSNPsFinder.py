import sys

LIST_PATH = "..\\data\\textFiles\\snpsToGene_DBP_LIST.txt"
MAPPING_PATH = "..\\data\\textFiles\\snpsToGene_DBP_ensemble.txt"
MISSING_PATH = "..\\data\\textFiles\\DBP_genes_missing_in_ensembl.txt"

if __name__ == "__main__":
    allSNPs = set()
    foundSNPs = set()

    with open(LIST_PATH, "r", encoding="utf-8") as listFile:
        lines = listFile.readlines()
        for snp in lines:
            allSNPs.add(snp.strip())
        
    with open(MAPPING_PATH, "r", encoding="utf-8") as mapFile:
        lines = mapFile.readlines()
        for line in lines:
            foundSNPs.add(line.strip().split(" ")[0].strip())

    missing = allSNPs.difference(foundSNPs)
    with open(MISSING_PATH, "w+", encoding="utf-8") as missingFile:
        for variant in missing:
            missingFile.write(variant + "\n")