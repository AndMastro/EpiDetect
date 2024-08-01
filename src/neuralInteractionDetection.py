import numpy as np
from itertools import combinations
from snpReader import snpRead
from snpReader import geneRead
from sklearn.metrics import roc_auc_score
import sys
from tqdm import tqdm
import os
import time

#INPUT_SIZE = 849  # 256*3  804 SBP, 1026 DBP, 849 PP
#DATASET_NAME = "SNPS_PP_AVG"
INPUT_SIZE = 849 #264*3 #2000*3 #256*3  # 256*3  804 SBP, 1026 DBP, 849 PP
DATASET_TYPE = "PP"
DATASET_NAME = "SNPS_PP_AVG"
if __name__ == "__main__":

    trueData = True  # variable stating if using true or generated data

    # if u wanna save on file
    save = False

    if sys.argv[1] == "save":
        save = True

    # specify degree of interaction

    interactionWay = int(sys.argv[2])

    layerWeights0 = None
    layerWeights1 = None
    layerWeights2 = None
    layerWeights3 = None
    sbpSNPs = None

    layerWeights0 = np.load(
        '../data/weights/BP/layer_0_weights_' + DATASET_NAME + '_numLayers2.npy', allow_pickle=True)
    layerWeights1 = np.load(
        '../data/weights/BP/layer_1_weights_' + DATASET_NAME + '_numLayers2.npy', allow_pickle=True)
    layerWeights2 = np.load(
        '../data/weights/BP/layer_2_weights_' + DATASET_NAME + '_numLayers2.npy', allow_pickle=True)
    layerWeights3 = np.load(
        '../data/weights/BP/layer_3_weights_' + DATASET_NAME + '_numLayers2.npy', allow_pickle=True)

    print("Loading SNPs database...")

    sbpSNPs = snpRead(
        "../data/snp_lists/allChrom_PP.bim")

    snpToGene = None

    if trueData == True:
        snpToGene = geneRead("../data/snp_gene_mappings/snpsToGene_" + DATASET_TYPE + "_ensembl_manual.txt",
                             "../data/snp_gene_mappings/notMapped_" + DATASET_TYPE + "_ensembl_closest.txt")

    '''print(layerWeights0)
    print(layerWeights1)
    print(layerWeights2)
    print(layerWeights3)'''

    dense1 = layerWeights0[0]
    dense2 = layerWeights1[0]
    output = layerWeights3[0]
    print("===========================")

    start = time.time()
    # print(output)

    # Define weights vectors according to NID paper
    # we need all the weights for all the input variables for any unit. So we need a matrix such thta we have a vector of 153 elems for any of the 64 hidden units

    print("Defining vectors for interaction detection. We need absolute values")
    # print(dense1[0])
    # print(len(dense1[0]))

    dense1Weights = []
    firstAppend = True

    for inputVec in dense1:
        for i in range(len(inputVec)):
            if firstAppend:
                dense1Weights.append([abs(inputVec[i])])
            else:
                dense1Weights[i].append(abs(inputVec[i]))
        firstAppend = False

    dense2Weights = []
    firstAppend = True

    for inputVec in dense2:
        for i in range(len(inputVec)):
            if firstAppend:
                dense2Weights.append([abs(inputVec[i])])
            else:
                dense2Weights[i].append(abs(inputVec[i]))
        firstAppend = False

    # print(dense1Weights[0])
    # print(len(dense1Weights[0]))

    # print(dense2Weights[0])
    # print(len(dense2Weights[0]))

    outputWeights = []
    for w in output:
        outputWeights.append((list(abs(w))))

    flatten = lambda l: [item for sublist in l for item in sublist]

    outputWeights = flatten(outputWeights)

    # print(len(outputWeights))

    print("Defining Aggregate Weight for any unit")

    # print(len(dense2))
    # print(len(output))

    aggregateWeightsD1 = np.matmul(
        np.array(outputWeights), np.array(dense2Weights))

    # print((aggregateWeightsD1))
    # print(len(aggregateWeightsD1))
    numInputs = list(np.arange(int(INPUT_SIZE)))

    potentialInteractions = list(combinations(numInputs, 2))

    # print((potentialInteractions))
    # print(len(potentialInteractions))
    print("There are " + str(len(potentialInteractions)) +
          " potential interactions.")
    
    print("Inizializing interaction dictionary.")
    interactions = {}
    
    for pair in potentialInteractions:
        interactions[pair] = []

    # print(interactions)

    # we use here aggregate weight and int avg at 1st layer (dens1weights)
    for intPair in tqdm(interactions):
        intStrength = 0
        for i in range(len(aggregateWeightsD1)):
            intStrength += aggregateWeightsD1[i]*(
                min(dense1Weights[i][intPair[0]], dense1Weights[i][intPair[1]]))
        interactions[intPair] = intStrength

    # print(interactions)

    # sort pairs according to interaction strenght
    intRank = sorted(((v, k) for k, v in interactions.items()), reverse=True)

    end = time.time()
    interactingSNPs = []

    if trueData:
        for pair in intRank:
            interactingSNPs.append((pair[0], (sbpSNPs[pair[1][0]//3], sbpSNPs[pair[1][1]//3]),
                                   (snpToGene[sbpSNPs[pair[1][0]//3]], snpToGene[sbpSNPs[pair[1][1]//3]])))  # add genes
    else:
        for pair in intRank:
                interactingSNPs.append(
                    (pair[0], (sbpSNPs[pair[1][0]//3], sbpSNPs[pair[1][1]//3])))
    # and now, analyze
        '''print(interactingSNPs[0:100])
        print(max(interactions.values()))
        print(min(interactions.values()))'''
    print("Time elapsed for NID: " + str(end-start))
    print("========================================================================")

    if save == True:
        counter = 0
        fileName = "../data/results/epistaticInteraction/NID/" + \
            str(interactionWay) + \
            "/thesis_11_10_23/epistaticInteractions_NID_" + DATASET_NAME + "_run{}.txt"
        while os.path.isfile(fileName.format(counter)):
            counter += 1
        fileName = fileName.format(counter)
        intFile = open(fileName, "w+")
        # intFile100 = open("../data/results/epistaticInteraction/268_SBP/prove/epistaticInteractions_268_200Units_NID_TOP_100_ONLYGENES.txt", "w+")

        # done in order to remove SNPs not related to genes (only kept for training since they intervene in SBP regulation but cannot raise epistatic phenomenon)
        if trueData:
            for pair in interactingSNPs:

                    snp1 = str(pair[1][0])
                    snp2 = str(pair[1][1])
                    gene1 = str(pair[2][0])
                    gene2 = str(pair[2][1])
                    score = str(pair[0])

                    intFile.write(
                        snp1 + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + score + "\n")
        else:
            for pair in interactingSNPs:
                snp1 = str(pair[1][0])
                snp2 = str(pair[1][1])
                score = str(pair[0])
                intFile.write(snp1 + "\t" + snp2 + "\t" + score + "\n")                    

        print("File with interactions saved.")

    sys.exit()