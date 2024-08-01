from collections import OrderedDict
import sys
if __name__ == "__main__":

    path = "../../data/results/epistaticInteraction/NID/BP/old_formats/epistaticInteractions_NID_SNPS_SBP_AVG_no_van_removed.txt"
    save_path = "../../data/results/epistaticInteraction/NID/BP/proper_format/epistaticInteractions_NID_SBP_SUM_ordered.txt"

    #simulated data
    '''interactions = {}
    with open(path, "r") as interFile:
        lines = interFile.readlines()
        for line in lines:
            line = line.strip().split("\t")
            snpPair = (line[0],line[1])
            score = line[2] 
            if snpPair not in interactions:
                interactions[snpPair] = float(score)
            else:
                interactions[snpPair] += float(score)
    
    intRank = sorted(((v, k) for k, v in interactions.items()), reverse=True)
    interactions = OrderedDict()
    for elem in intRank:
        interactions[elem[1]] = elem[0]

    with open(save_path, "w+") as saveFile:
        for pair in interactions:
            snp1 = pair[0]
            snp2 = pair[1]
            score = interactions[pair]
            if snp1 != snp2:
                saveFile.write(snp1 + "\t" + snp2 + "\t" + str(score) + "\n")'''

    #real data

    interactions = {}
    with open(path, "r") as interFile:
        lines = interFile.readlines()
        for line in lines:
            line = line.strip().split("\t")
            snpPair = (line[0],line[1], line[2], line[3])
            score = line[4] 
            if snpPair not in interactions:
                interactions[snpPair] = float(score)
            else:
                interactions[snpPair] += float(score)
    
    interactions = dict(sorted(interactions.items(), reverse = True, key=lambda item: item[1]))
    
    
    with open(save_path, "w+") as saveFile:
        for pair in interactions:
            snp1 = pair[0]
            snp2 = pair[1]
            gene1 = pair[2]
            gene2 = pair[3]
            score = interactions[pair]
            if snp1 != snp2:
                saveFile.write(snp1 + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + str(score) + "\n")
        