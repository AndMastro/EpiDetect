

def snpRead(path):
    fp = open(path, 'r')
    numSNPs = 0
    lines = fp.readlines()
    sbpSNPs = {}
    for line in lines:
        sbpSNPs[numSNPs] = line.split("\t")[1]
        numSNPs+=1
        
    return sbpSNPs

def geneRead(path):
    fp = open(path, 'r')
    lines = fp.readlines()
    snpToGene = {}

    for line in lines:
        line = line.strip().split(" ")
        snpID = line[0].replace(" ","")
        gene = line[-1].replace(" ","")
        snpToGene[snpID] = gene
        
        
    return snpToGene

if __name__ == "__main__":

    snps = snpRead("F:\\Repositories\\BioNet\\data\\binFiles\\264SNPs\\allChromImp.bim")
    print(snps)

    snpToGene = geneRead("F:\\Repositories\\BioNet\\data\\textFiles\\snpsToGene.txt")
    print(snpToGene["rs7187540"])