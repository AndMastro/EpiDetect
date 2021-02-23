from snpReader import geneRead, snpRead
import sys
from tqdm import tqdm 

#this file converst the ME output of BOOST into a version in which we can see the SNP rsID and the gene name (and not the index)
if __name__ == "__main__":
    BP = "PP"
    PATH = "../data/results/epistaticInteraction/MDR/BP/old_format/marginal_effect_MDR_" + BP + "_SNPs_only.txt"
    SAVE_PATH = "../data/results/epistaticInteraction/MDR/BP/proper_format/marginal_effect_MDR_" + BP + ".txt"

    snpDict = sbpSNPs = snpRead("../data/snp_lists/allChrom_" + BP +".bim")
    
    snpToGene = geneRead("../data/snp_gene_mappings/snpsToGene_" + BP +"_ensembl_manual.txt",
                             "../data/snp_gene_mappings/notMapped_" + BP +"_ensembl_closest.txt")


    save_file = open(SAVE_PATH, "w+")

    marginal_effects = {}
    with open(PATH, "r") as interFile:
        lines = interFile.readlines()
        for line in lines:
            #print(line)
            #sys.exit()
            line = line.strip().split("\t")
            snp = line[0].replace(" ", "")
            score = line[1].replace(" ", "") 
            
            
            #marginal_effects[snp] = abs(float(score))
            marginal_effects[snp] = float(score)
    
    # meRank = sorted(((v, k) for k, v in marginal_effects.items()), reverse=True)
    meRank = sorted(((v, k) for k, v in marginal_effects.items()), reverse=True) #biased toward high marginal effect with negative effect. Makes sense tbh, very strong negative ME.
    marginal_effects = {}
    for elem in meRank:
        marginal_effects[elem[1]] = elem[0]

    #print(marginal_effects)
    
    with open(SAVE_PATH, "w+") as saveFile:
        for snp in marginal_effects:
            #snp1 = snpDict[int(snp)]
            snp1 = snp
            #snp2 = snpDict[int(pair[1])]
            score = marginal_effects[snp]

            gene1 = snpToGene[snp1]
            #gene2 = snpToGene[snp2]
            # print(snp1)
            # print(gene1)
            #sys.exit()
            save_file.write(snp1  + "\t" + gene1 + "\t" + str(score) + "\n")
