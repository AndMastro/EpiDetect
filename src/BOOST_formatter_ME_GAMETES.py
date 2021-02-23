from snpReader import geneRead, snpRead
import sys
from tqdm import tqdm 

#this file converst the ME output of BOOST into a version in which we can see the SNP rsID and the gene name (and not the index)
if __name__ == "__main__":
    DS = "epistatic"
    PATH = "../data/results/epistaticInteraction/BOOST/GAMETES/" + DS + "/BOOST_marginal_effects_epistatic_MAF01_eta01_theta1.txt"
    SAVE_PATH = "../data/results/epistaticInteraction/BOOST/GAMETES/" + DS + "/BOOST_marginal_effects_epistatic_MAF01_eta01_theta1_sorted.txt"

    snpDict = sbpSNPs = snpRead("../data/snp_lists/snpList256.bim")
    

    save_file = open(SAVE_PATH, "w+")

    marginal_effects = {}
    with open(PATH, "r") as interFile:
        lines = interFile.readlines()
        for line in lines:
            #print(line)
            #sys.exit()
            line = line.strip().split(" ")
            snp = line[0].replace(" ", "")
            score = line[1].replace(" ", "") 
            
            
            #marginal_effects[snp] = abs(float(score))
            marginal_effects[snp] = float(score)
    
    # meRank = sorted(((v, k) for k, v in marginal_effects.items()), reverse=True)
    meRank = sorted(((v, k) for k, v in marginal_effects.items()), reverse=True)
    marginal_effects = {}
    for elem in meRank:
        marginal_effects[elem[1]] = elem[0]

    #print(marginal_effects)
    
    with open(SAVE_PATH, "w+") as saveFile:
        for snp in marginal_effects:
            snp1 = snpDict[int(snp)]
            #snp2 = snpDict[int(pair[1])]
            score = marginal_effects[snp]

            #gene2 = snpToGene[snp2]
            # print(snp1)
            # print(gene1)
            #sys.exit()
            save_file.write(snp1 + "\t" + str(score) + "\n")
