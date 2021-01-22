import sys

METHOD  = "BOOST"
TRAIT = "SBP"
PATH = "../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ METHOD +"_"+TRAIT+ ".txt"

SNP_genes_dict = {}
with open (PATH, "r") as interactionFile:
    lines = interactionFile.readlines()[:100]
    for line in lines:
        row = line.strip().split("\t")
        snp1 = row[0]
        snp2 = row[1]
        gene1 = row[2]
        gene2 = row[3]
        
        if (snp1, gene1) not in SNP_genes_dict:
            SNP_genes_dict[(snp1, gene1)] = 1
        else:
            SNP_genes_dict[(snp1, gene1)] += 1

        if (snp2, gene2) not in SNP_genes_dict:
            SNP_genes_dict[(snp2, gene2)] = 1
        else:
            SNP_genes_dict[(snp2, gene2)] += 1

top = sorted(((v, k) for k, v in SNP_genes_dict.items()), reverse=True)
print(top[:5])
print("Num SNPs in top1000 interactions: " +  str(len(top)))
