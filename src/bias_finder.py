import sys


top_n_list = [100, 500, 1000]

for top_n in top_n_list:

    TRAIT = "PP"
    #BOOST

    METHOD  = "BOOST"

    PATH = "../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ METHOD +"_"+TRAIT+ ".txt" 

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
    #print("Num SNPs in top1000 interactions: " +  str(len(top)))

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
    #print("BOOST:")
    #print(snps)
    BOOST_top = [elem[0] for elem in snps]



    #MDR

    METHOD  = "MDR"

    PATH = "../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/epistaticInteractions_"+ TRAIT +"_"+METHOD+ "_top1000.txt"

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
    #print("Num SNPs in top1000 interactions: " +  str(len(top)))


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
    MDR_top = [elem[0] for elem in snps]

    METHOD  = "NID"

    PATH = "../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/thesis_11_10_2023/epistaticInteractions_"+ METHOD +"_"+TRAIT+ "_SUM_ordered.txt"


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
    #print("Num SNPs in top1000 interactions: " +  str(len(top)))

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
    NID_top = [elem[0] for elem in snps]

    #EpiCID thesis contatins already no van removed

    METHOD  = "EpiCiD"

    PATH = "../data/results/epistaticInteraction/"+ METHOD +"/BP/proper_format/thesis_11_10_2023/epistaticInteractions_"+ METHOD +"_"+TRAIT+ ".txt"

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
    #print("Num SNPs in top1000 interactions: " +  str(len(top)))


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

    #SBP
    #BOOST VS EPICID VS MDR TOP 100

    # plt.title('Most interacting SNPS - Top 100 interactions')
    # p1 = plt.bar([1,2,3,4], [100, 78, 35,0])
    # p2 = plt.bar([1,2,3,4], [0,19, 27,0], bottom = [100,78,35,0])
    # p3 = plt.bar([1,2,3,4], [0,3,18,0], bottom = [0,78 +19, 27+35,0])
    # p4 = plt.bar([1,2,3,4], [0,0,11,0], bottom = [0,78 +19,27+35+18,0])
    # p5 = plt.bar([1,2,3,4], [0,0,9,0], bottom = [0,78 +19, 27+35+18+11,0])
    # plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Top 1 SNP', 'Top 2 SNP', 'Top 3 SNP', 'Top 4 SNP', 'Top 5 SNP'))
    # plt.xticks([1,2,3], ('BOOST', 'MDR', 'EpiCID'))
    # plt.ylabel('Number of interactions')
    # plt.xlabel('Method')
    # plt.show()


    #BOOST VS EPICID VS MDR TOP 500

    # plt.title('Most interacting SNPS - Top 500 interactions')
    # p1 = plt.bar([1,2,3,4], [216,267,96,0])
    # p2 = plt.bar([1,2,3,4], [182,186,85,0], bottom = [216,267,96,0])
    # p3 = plt.bar([1,2,3,4], [56,32,59,0], bottom = [216+182,267+186,96+85,0])
    # p4 = plt.bar([1,2,3,4], [46,15,41,0], bottom = [216+182+56,267+186+32,96+85+59,0])
    # p5 = plt.bar([1,2,3,4], [0,0,41,0], bottom = [216+182+56,267+186+32,96+85+59+41,0])
    # p6 = plt.bar([1,2,3,4], [0,0,33,0], bottom = [216+182+56,267+186+32,96+85+59+41+41,0])
    # p7 = plt.bar([1,2,3,4], [0,0,31,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33,0])
    # p8 = plt.bar([1,2,3,4], [0,0,28,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31,0])
    # p9 = plt.bar([1,2,3, 4], [0,0,21,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31+28,0])
    # p10 = plt.bar([1,2,3,4], [0,0,19,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31+28+21,0])
    # p11 = plt.bar([1,2,3,4], [0,0,19,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31+28+21+19,0])
    # p12 = plt.bar([1,2,3,4], [0,0,18,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31+28+21+19+19,0])
    # p12 = plt.bar([1,2,3,4], [0,0,9,0], bottom = [216+182+56,267+186+32,96+85+59+41+41+33+31+28+21+19+19+18,0])
    # plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0], p9[0], p10[0], p11[0], p12[0]), ('Top 1 SNP', 'Top 2 SNP', 'Top 3 SNP', 'Top 4 SNP', 'Top 5 SNP', 'Top 6 SNP', 'Top 7 SNP', 'Top 8 SNP', 'Top 9 SNP', 'Top 10 SNP','Top 11 SNP','Top 12 SNP'))
    # plt.xticks([1,2,3], ('BOOST', 'MDR', 'EpiCID'))
    # plt.ylabel('Number of interactions')
    # plt.xlabel('Method')
    # plt.show()


    #BOOST vs EpiCID vs MDR top 1000

    # plt.title('Most interacting SNPS - Top 1000 interactions')
    # p1 = plt.bar([1,2,3,4], [230,267,139,0])
    # p2 = plt.bar([1,2,3,4], [182,267,124,0], bottom = [230,267,139,0])
    # p3 = plt.bar([1,2,3,4], [160,239,89,0], bottom = [230+182,267+267,139+124,0])
    # p4 = plt.bar([1,2,3,4], [132,156,74,0], bottom = [230+182+160,267+267+239,139+124+89,0])
    # p5 = plt.bar([1,2,3,4], [91,45,67,0], bottom = [230+182+160+132,267+267+239+156,139+124+89+74,0])
    # p6 = plt.bar([1,2,3,4], [59,26,63,0], bottom = [230+182+160+132+91,267+267+239+156+45,139+124+89+74+67,0])
    # p7 = plt.bar([1,2,3,4], [55,0,52,0], bottom = [230+182+160+132+91+59,267+267+239+156+45,139+124+89+74+67+63,0])
    # p8 = plt.bar([1,2,3,4], [51,0,48,0], bottom = [230+182+160+132+91+59+55,267+267+239+156+45,139+124+89+74+67+63+52,0])
    # p9 = plt.bar([1,2,3,4], [40,0,41,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48,0])
    # p10 = plt.bar([1,2,3,4], [0,0,38,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41,0])
    # p11 = plt.bar([1,2,3,4], [0,0,36,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38,0])
    # p12 = plt.bar([1,2,3,4], [0,0,33,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36,0])
    # p13 = plt.bar([1,2,3,4], [0,0,33,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33,0])
    # p14 = plt.bar([1,2,3,4], [0,0,32,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33,0])
    # p15 = plt.bar([1,2,3,4], [0,0,32,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33+32,0])
    # p16 = plt.bar([1,2,3,4], [0,0,30,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33+32+32,0])
    # p17 = plt.bar([1,2,3,4], [0,0,27,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33+32+32+30,0])
    # p18 = plt.bar([1,2,3,4], [0,0,26,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33+32+32+30+27,0])
    # p19 = plt.bar([1,2,3,4], [0,0,16,0], bottom = [230+182+160+132+91+59+55+51,267+267+239+156+45,139+124+89+74+67+63+52+48+41+38+36+33+33+32+32+30+27+26,0])
    # plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0], p9[0], p10[0], p11[0], p12[0], p13[0],p14[0],p15[0],p16[0],
    #                 p17[0],p18[0],p19[0]), ('Top 1 SNP', 'Top 2 SNP', 'Top 3 SNP', 'Top 4 SNP', 'Top 5 SNP', 'Top 6 SNP', 
    #                 'Top 7 SNP', 'Top 8 SNP', 'Top 9 SNP', 'Top 10 SNP','Top 11 SNP','Top 12 SNP'
    #                             ,'Top 13 SNP','Top 14 SNP','Top 15 SNP','Top 16 SNP','Top 17 SNP','Top 18 SNP','Top 19 SNP'))
    # plt.xticks([1,2,3], ('BOOST', 'MDR', 'EpiCID'))
    # plt.ylabel('Number of interactions')
    # plt.xlabel('Method')
    # plt.show()

    #DBP

    #BOOST VS MDR VS NID vs EpiCID TOP 100

    # plt.title('Most interacting SNPs - Top 100 interactions')
    # p1 = plt.bar([1,2,3,4,5], [100, 45, 28, 25, 0])
    # p2 = plt.bar([1,2,3,4,5], [0,27,14,19,0], bottom = [100, 45, 28, 25, 0])
    # p3 = plt.bar([1,2,3,4,5], [0,18,13,15,0], bottom = [100,45 +27, 28+14, 25+19,0])
    # p4 = plt.bar([1,2,3,4,5], [0,10,11,15,0], bottom = [100,45 +27+18, 28+14+13, 25+19+15,0])
    # p5 = plt.bar([1,2,3,4,5], [0,0,11,13,0], bottom = [100,45 +27+18+10, 28+14+13+11, 25+19+15+15,0])
    # p6 = plt.bar([1,2,3,4,5], [0,0,11,12,0], bottom = [100,45 +27+18+10, 28+14+13+11+11, 25+19+15+15+13,0])
    # p7 = plt.bar([1,2,3,4,5], [0,0,10,1,0], bottom = [100,45 +27+18+10, 28+14+13+11+11+11, 25+19+15+15+13+12,0])
    # p8 = plt.bar([1,2,3,4,5], [0,0,2,0,0], bottom = [100,45 +27+18+10, 28+14+13+11+11+11+10, 25+19+15+15+13+12,0])

    # plt.legend((p1[0], p2[0], p3[0], p4[0],p5[0],p6[0],p7[0],p8[0]), ('Top 1 SNP', 'Top 2 SNP', 'Top 3 SNP', 'Top 4 SNP', 'Top 5 SNP', 'Top 6 SNP', 'Top 7 SNP', 'Top 8 SNP'))
    # plt.xticks([1,2,3,4], ('BOOST', 'MDR', 'NID', 'EpiCID'))
    # plt.ylabel('Number of interactions')
    # plt.xlabel('Method')
    # plt.show()

    #BOOST VS MDR VS NID vs EpiCID TOP 500

    # plt.title('Most interacting SNPs - Top 500 interactions')
    # p1 = plt.bar([1,2,3,4,5], [246, 160, 76, 84, 0])
    # p2 = plt.bar([1,2,3,4,5], [181,114,45,51,0], bottom = [246, 160, 76, 84, 0])
    # p3 = plt.bar([1,2,3,4,5], [73,77,44,51,0], bottom = [246+181, 160+114, 76+45, 84+51, 0])
    # p4 = plt.bar([1,2,3,4,5], [0,72,43,47,0], bottom = [246+181+73, 160+114+77, 76+45+44, 84+51+51, 0])
    # p5 = plt.bar([1,2,3,4,5], [0,40,43,46,0], bottom = [246+181+73, 160+114+77+72, 76+45+44+43, 84+51+51+47, 0])
    # p6 = plt.bar([1,2,3,4,5], [0,27,42,45,0], bottom = [246+181+73, 160+114+77+72+40, 76+45+44+43+43, 84+51+51+47+46, 0])
    # p7 = plt.bar([1,2,3,4,5], [0,10,35,35,0], bottom = [246+181+73, 160+114+77+72+40+27, 76+45+44+43+43+42, 84+51+51+47+46+45, 0])
    # p8 = plt.bar([1,2,3,4,5], [0,0,35,31,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35, 84+51+51+47+46+45+35, 0])
    # p9 = plt.bar([1,2,3,4,5], [0,0,34,21,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35, 84+51+51+47+46+45+35+31, 0])
    # p10 = plt.bar([1,2,3,4,5], [0,0,29,21,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35+34, 84+51+51+47+46+45+35+31+21, 0])
    # p11 = plt.bar([1,2,3,4,5], [0,0,27,20,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35+34+29, 84+51+51+47+46+45+35+31+21+21, 0])
    # p12 = plt.bar([1,2,3,4,5], [0,0,26,19,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35+34+29+27, 84+51+51+47+46+45+35+31+21+21+20, 0])
    # p13 = plt.bar([1,2,3,4,5], [0,0,21,19,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35+34+29+27+26, 84+51+51+47+46+45+35+31+21+21+20+19, 0])
    # p14 = plt.bar([1,2,3,4,5], [0,0,0,10,0], bottom = [246+181+73, 160+114+77+72+40+27+10, 76+45+44+43+43+42+35+35+34+29+27+26+21, 84+51+51+47+46+45+35+31+21+21+20+19+19, 0])

    # sub_plots= [p1, p2 ,p3 ,p4 ,p5 ,p6 ,p7 ,p8 ,p9 ,p10,p11,p12,p13,p14]
    # x_legend = [p[0] for p in sub_plots]
    # name_legend = []
    # for i in range(1, 14 + 1):
    #     name_legend.append("Top " + str(i) + " SNP")

    # plt.legend(x_legend, name_legend)
    # plt.xticks([1,2,3,4], ('BOOST', 'MDR', 'NID', 'EpiCID'))
    # plt.ylabel('Number of interactions')
    # plt.xlabel('Method')
    # plt.show()

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

    plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree_thesis_28_10_23.png", dpi = 300, bbox_inches='tight')
    plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree_thesis_28_10_23.svg", dpi = 300, bbox_inches='tight')
    plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree_thesis_28_10_23.eps", dpi = 300, bbox_inches='tight')
    plt.savefig("../data/results/epistaticInteraction/comparison_plots/" + TRAIT + "/epistaticInteractions_top" + str(top_n) + "_highest_degree_thesis_28_10_23.pdf", dpi = 300, bbox_inches='tight')
