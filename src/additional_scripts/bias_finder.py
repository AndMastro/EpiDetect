import sys


top_n_list = [100, 500, 1000]

for top_n in top_n_list:

    TRAIT = "PP"
    #BOOST

    METHOD  = "BOOST"

    PATH = "../../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ METHOD +"_"+TRAIT+ ".txt" 

    SNP_genes_dict = {}
    with open (PATH, "r") as interactionFile:
        lines = interactionFile.readlines()[:top_n]
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
    print("BOOST:")
    print(top[:5])
    

    #find most interacting SNPs among first top_n
    snps = []
    count = 0
    for elem in top:
        count += elem[0]
        snps.append(elem)
        if count >= top_n:
            break
    
    sum_occurences = 0

    for elem in snps:
        sum_occurences += elem[0]

   
    snps[-1] = (snps[-1][0] - (sum_occurences  - top_n), (snps[-1][1][0], snps[-1][1][1]))
    
    BOOST_top = [elem[0] for elem in snps]



    #MDR

    METHOD  = "MDR"

    PATH = "../../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ TRAIT +"_"+METHOD+ "_top1000.txt"

    SNP_genes_dict = {}
    with open (PATH, "r") as interactionFile:
        lines = interactionFile.readlines()[:top_n]
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
    print("MDR:")
    print(top[:5])
    


    #find most interacting SNPs among first top_n
    snps = []
    count = 0
    for elem in top:
        count += elem[0]
        snps.append(elem)
        if count >= top_n:
            break
    
    sum_occurences = 0

    for elem in snps:
        sum_occurences += elem[0]

    
    snps[-1] = (snps[-1][0] - (sum_occurences  - top_n), (snps[-1][1][0], snps[-1][1][1]))
    
    MDR_top = [elem[0] for elem in snps]

    METHOD  = "NID"

    PATH = "../../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ METHOD +"_"+TRAIT+ "_SUM_ordered.txt"


    #NID
    SNP_genes_dict = {}
    with open (PATH, "r") as interactionFile:
        lines = interactionFile.readlines()[:top_n]
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
    print("NID:")
    print(top[:5])
    

    #find most interacting SNPs among first top_n
    snps = []
    count = 0
    for elem in top:
        count += elem[0]
        snps.append(elem)
        if count >= top_n:
            break
    
    sum_occurences = 0

    for elem in snps:
        sum_occurences += elem[0]

    
    snps[-1] = (snps[-1][0] - (sum_occurences  - top_n), (snps[-1][1][0], snps[-1][1][1]))

    
    NID_top = [elem[0] for elem in snps]

    

    METHOD  = "EpiCiD"

    PATH = "../../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ METHOD +"_"+TRAIT+ ".txt"

    SNP_genes_dict = {}
    with open (PATH, "r") as interactionFile:
        lines = interactionFile.readlines()[:top_n]
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
    print("EpiCID")
    print(top[:5])
    


    #find most interacting SNPs among first top_n
    snps = []
    count = 0
    for elem in top:
        count += elem[0]
        snps.append(elem)
        if count >= top_n:
            break
    #print(snps)
    sum_occurences = 0

    for elem in snps:
        sum_occurences += elem[0]

    #print(sum_occurences)
    snps[-1] = (snps[-1][0] - (sum_occurences  - top_n), (snps[-1][1][0], snps[-1][1][1]))

    #print(snps)
    EpiCID_top = [elem[0] for elem in snps]

    import numpy as np
    import matplotlib.pyplot as plt

    

    #BOOST VS MDR VS NID vs EpiCID GENERAL

    max_snps = max([len(BOOST_top), len(MDR_top), len(NID_top), len(EpiCID_top)])

    BOOST_top = BOOST_top + [0]*(max_snps - len(BOOST_top))
    MDR_top = MDR_top + [0]*(max_snps - len(MDR_top))
    NID_top = NID_top + [0]*(max_snps - len(NID_top))
    EpiCID_top = EpiCID_top + [0]*(max_snps - len(EpiCID_top))

    plt.rcParams.update({'font.size': 23})
    plt.figure(figsize= (14,12))
    # plt.tight_layout()
    plt.grid(visible=True, linestyle='-', axis="y", linewidth=0.5, alpha=0.8)
    # plt.title('Highest-degree SNPs - Top-' + str(top_n) +' interactions')
    sub_plots = []

    # p = plt.bar([1,2,3,4,5], [BOOST_top[0], MDR_top[0], NID_top[0],EpiCID_top[0],0])
    p = plt.bar([1,2,3,4], [BOOST_top[0], MDR_top[0], NID_top[0],EpiCID_top[0]])
    sub_plots.append(p)

    for i in range(1, max_snps):
        # y_bar = [BOOST_top[i], MDR_top[i], NID_top[i],EpiCID_top[i],0]
        y_bar = [BOOST_top[i], MDR_top[i], NID_top[i],EpiCID_top[i]]
        # bottom_vals = [sum(BOOST_top[0:i]), sum(MDR_top[0:i]),sum(NID_top[0:i]),sum(EpiCID_top[0:i]),0]
        bottom_vals = [sum(BOOST_top[0:i]), sum(MDR_top[0:i]),sum(NID_top[0:i]),sum(EpiCID_top[0:i])]
        # p = plt.bar([1,2,3,4,5], y_bar, bottom=bottom_vals)
        p = plt.bar([1,2,3,4], y_bar, bottom=bottom_vals)
        sub_plots.append(p)


    x_legend = [p[0] for p in sub_plots]
    name_legend = []
    for i in range(1, max_snps + 1):
        # name_legend.append("Top-" + str(i) + " SNP")
        name_legend.append(str(i))

    plt.legend(x_legend, name_legend, title="SNP rank",bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.xticks([1,2,3,4], ('BOOST', 'MDR', 'NID', 'EpiCID'))
    plt.ylabel('Number of interactions')
    plt.xlabel('Method')
    #plt.ylim(top=top_n+0.05*top_n) 
    plt.ylim(top=top_n) 
    # plt.show()

    plt.tight_layout()

    plt.savefig("../../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree.png", dpi = 300, bbox_inches='tight')
    plt.savefig("../../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree.svg", dpi = 300, bbox_inches='tight')
    plt.savefig("../../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree.eps", dpi = 300, bbox_inches='tight')
    plt.savefig("../../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree.pdf", dpi = 300, bbox_inches='tight')
