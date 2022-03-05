import csv
import sys

gene_list = []
gene_to_snp = {}
gene_scores = {}
pair_scores = {}
with open("marginal_effect_SBP_MDR_no_van_removed.txt", "r") as csv_file:
    reader = csv.reader(csv_file, delimiter = "\t")
    for line in reader:
        snp = line[0]
        gene = line[1]
        score = line[2]
        gene_list.append(gene)
        gene_to_snp[gene] = snp
        gene_scores[gene] = score

for gene1 in gene_list:
    for gene2 in gene_list:
        if gene1 != gene2:
            if (gene1, gene2) not in pair_scores and (gene2, gene1) not in pair_scores: #we avoid to add pairs twice
                pair_scores[(gene2, gene1)] = round(float(gene_scores[gene1]) + float(gene_scores[gene2]), 6)

pair_rankings = sorted(pair_scores, key=pair_scores.get, reverse=True) 

with open("additive_ranking_SBP_no_van_removed.txt", "w+") as csv_file:
    # writer = csv.writer(csv_file)
    for pair in pair_rankings:
        gene1 = pair[0]
        gene2 = pair[1]
        snp1 = gene_to_snp[gene1]
        snp2 = gene_to_snp[gene2]
        additive_score = str(pair_scores[pair])
        csv_file.write(snp1 + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + additive_score + "\n")