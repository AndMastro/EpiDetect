from snpReader import geneRead

if __name__ == "__main__":
    snpToGene = geneRead("..\\data\\snp_gene_mappings\\snpsToGene_SBP_ensembl_manual.txt", "..\\data\\snp_gene_mappings\\notMapped_SBP_ensembl_closest.txt")