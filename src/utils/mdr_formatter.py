from snpReader import geneRead
import sys
from tqdm import tqdm 

if __name__ == "__main__":

    realData = sys.argv[1]

    BP = "SBP"
    PATH = "../data/results/epistaticInteraction/MDR/BP/" + BP + "_AVG_analysis_MDR_no_van_removed.txt"
    SAVE_PATH = "../data/results/epistaticInteraction/MDR/"+ BP + "_AVG_MDR_genes_rank_top1000_no_van_removed.txt"

    if realData == "true":
    

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

    elif realData == "false":
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

    else:
        print("FATAL ERROR: Command not recognized. Aborting.")