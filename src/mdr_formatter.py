from snpReader import geneRead
import sys
from tqdm import tqdm 

if __name__ == "__main__":
    BP = "PP"
    PATH = "../data/results/epistaticInteraction/MDR/" + BP + "_AVG_MDR_rank.txt"
    SAVE_PATH = "../data/results/epistaticInteraction/MDR/"+ BP + " _AVG_MDR_genes_rank.txt"

    snpToGene = geneRead("../data/snp_gene_mappings/snpsToGene_" + BP +"_ensembl_manual.txt",
                             "../data/snp_gene_mappings/notMapped_" + BP +"_ensembl_closest.txt")


    save_file = open(SAVE_PATH, "w+")

    with open(PATH, "r") as rankFile:
        lines = rankFile.readlines()
        for line in tqdm(lines):
            row = line.strip().split(" ")
            snp1 = row[0]
            snp2 = row[1]
            score = row[2]
            
            gene1 = snpToGene[snp1]
            gene2 = snpToGene[snp2]

            save_file.write(snp1  + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + score + "\n")

    save_file.close()