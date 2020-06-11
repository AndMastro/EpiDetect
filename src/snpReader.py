

def snpRead(path):
    fp = open(path, 'r')
    numSNPs = 0
    lines = fp.readlines()
    sbpSNPs = {}
    for line in lines:
        sbpSNPs[numSNPs] = line.split("\t")[1]
        numSNPs+=1
        
    return sbpSNPs

def snpRead_chr_pos(path):
    fp = open(path, 'r')
    numSNPs = 0
    lines = fp.readlines()
    sbpSNPs = {}
    for line in lines:
        chrom = line.split("\t")[0]
        snp = line.split("\t")[1]
        pos = line.split("\t")[3]
        sbpSNPs[numSNPs] = [snp, chrom, pos]
        numSNPs+=1
        
    return sbpSNPs

def geneRead(path, path_missing):
    fp = open(path, 'r')
    lines = fp.readlines()
    snpToGene = {}

    for line in lines:
        line = line.strip().split(" ")
        snpID = line[0].replace(" ","")
        gene = line[-1].replace(" ","")
        snpToGene[snpID] = gene

    fp.close()
        
    fpm = open(path_missing, 'r')
    lines = fpm.readlines()
    
    missing_manual_mapping = {}

    for line in lines:
        line = line.strip().split(" ")
        snpID = line[0].replace(" ","")
        gene = line[-1].replace(" ","")
        missing_manual_mapping[snpID] = gene
    
    #putting manually annotate gene in lieu of None

    for snp in snpToGene:
        if snpToGene[snp] == "None":
            snpToGene[snp] = missing_manual_mapping[snp]

    fpm.close()

    return snpToGene

if __name__ == "__main__":

    snps = snpRead("F:\\Repositories\\BioNet\\data\\binFiles\\264SNPs\\allChromImp.bim")
    print(snps)

    snpToGene = geneRead("F:\\Repositories\\BioNet\\data\\textFiles\\snpsToGene.txt")
    print(snpToGene["rs7187540"])