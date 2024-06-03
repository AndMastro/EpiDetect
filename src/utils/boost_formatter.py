from snpReader import geneRead, snpRead
import sys
from tqdm import tqdm 

if __name__ == "__main__":
    BP = "SBP"
    PATH = "../data/results/epistaticInteraction/BOOST/BP/old_format/InteractionRecords_" + BP + "_no_van_removed.txt"
    SAVE_PATH = "../data/results/epistaticInteraction/BOOST/BP/proper_format/epistaticInteractions_BOOST_pairwise_" + BP + "_no_van_removed.txt"

    snpDict = sbpSNPs = snpRead("../data/snp_lists/allChrom_" + BP +"_no_van_removed.bim")
    
    snpToGene = geneRead("../data/snp_gene_mappings/snpsToGene_" + BP +"_ensembl_manual.txt",
                             "../data/snp_gene_mappings/notMapped_" + BP +"_ensembl_closest.txt")


    save_file = open(SAVE_PATH, "w+")

    interactions = {}
    with open(PATH, "r") as interFile:
        lines = interFile.readlines()
        for line in lines:
            line = line.strip().split("\t")
            snpPair = (line[1].replace(" ", ""),line[2].replace(" ", ""))
            score = line[5].replace(" ", "") 
            
            interactions[snpPair] = float(score)

    intRank = sorted(((v, k) for k, v in interactions.items()), reverse=True)
    interactions = {}
    for elem in intRank:
        interactions[elem[1]] = elem[0]

    with open(SAVE_PATH, "w+") as saveFile:
        for pair in interactions:
            snp1 = snpDict[int(pair[0])]
            snp2 = snpDict[int(pair[1])]
            score = interactions[pair]

            gene1 = snpToGene[snp1]
            gene2 = snpToGene[snp2]

            save_file.write(snp1  + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + str(score) + "\n")
