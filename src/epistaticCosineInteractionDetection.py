import numpy as np
import os
from itertools import combinations
from snpReader import snpRead
from snpReader import geneRead
from sklearn.metrics import pairwise
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

INPUT_SIZE =  849 #264*3 #2000*3 #256*3  # 264*3  804 SBP, 1026 DBP, 849 PP
DATASET_NAME = "SNPS_PP_AVG"
DATASET_TYPE =  "PP"
# DATASET_TYPE = "marginal_effect" 
# DATASET_NAME = "epistatic_plus_ME_risk_model_MAF01_eta01_theta1_lambda03_1_2_EDM-1_01"  #used only with GAMETES data
if __name__ == "__main__":

    trueData = True  # variable stating if using true or generated data by GAMETES

    # if you want save on file
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
        "../data/snp_lists/allChrom_PP.bim") #allChrom_SBP #allChrom_PP #allChrom_DBP

    snpToGene = None

    if trueData == True:
        snpToGene = geneRead("../data/snp_gene_mappings/snpsToGene_" + DATASET_TYPE +"_ensembl_manual.txt",
                             "../data/snp_gene_mappings/notMapped_" + DATASET_TYPE +"_ensembl_closest.txt")

    dense1 = layerWeights0[0]
    dense2 = layerWeights1[0]
    output = layerWeights3[0]

    print("Done!")
    print("===========================")
    
    start = time.time()

    # Define weight vectors 
    print("Defining vectors for interaction detection. We need absolute values...")

    dense1Weights = []
    firstAppend = True

    # for any input unit we have a vector with the weights to all the fist layer hidden units
    for vector in dense1:
        dense1Weights.append(np.array(list(abs(vector))))

    dense2Weights = []
    firstAppend = True

    for inputVec in dense2:
        for i in range(len(inputVec)):
            if firstAppend:
                dense2Weights.append([abs(inputVec[i])])
            else:
                dense2Weights[i].append(abs(inputVec[i]))
        firstAppend = False

    outputWeights = []
    for w in output:
        outputWeights.append((list(abs(w))))

    def flatten(l): return [item for sublist in l for item in sublist]

    outputWeights = flatten(outputWeights)

    print("Defining Aggregated Weight for any unit...")

    aggregateWeightsD1 = np.matmul(
        np.array(outputWeights), np.array(dense2Weights))

    print("Done.")

    # defining aggregated neural feature vector for any SNP (aggregation of 3 input features as follows)
    numSNPs = 0
    snpsDict = {}
    scalingFactor11 = 0.5
    scalingFactor12 = 1
    scalingFactor22 = 1 
    for i in tqdm(range(0, len(dense1Weights), 3)):
        # uncomment accroding to if you want NFV sum OR concatenation
        neuralFeatureVectorProxy1 = np.multiply(
            dense1Weights[i], aggregateWeightsD1)
        neuralFeatureVectorProxy2 = np.multiply(
            dense1Weights[i+1], aggregateWeightsD1)
        neuralFeatureVectorProxy3 = np.multiply(
            dense1Weights[i+2], aggregateWeightsD1)

        # uncomment those two lines to have NFV sum
        neuralFeatureVector = np.add(
            neuralFeatureVectorProxy1, neuralFeatureVectorProxy2)
        neuralFeatureVector = np.add(
            neuralFeatureVector, neuralFeatureVectorProxy3) #*scalingFactor22

        # uncomment this line to have NFV concatenation
        #neuralFeatureVector = np.concatenate((neuralFeatureVectorProxy1, neuralFeatureVectorProxy2, neuralFeatureVectorProxy3))

        # uncomment this to have max
        #neuralFeatureVector = []
        '''for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(max(neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]))
        minVector = []'''

        # for min
        '''neuralFeatureVector = []
        for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(min(neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]))
        minVector = []'''

        # for avg
        '''neuralFeatureVector = []
        for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(np.mean(np.array((neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]))))
        minVector = []'''

        minVector = []
        concat = False

        if concat:
            minVector = np.concatenate(
                (dense1Weights[i], dense1Weights[i+1], dense1Weights[i+2]))
        else:
            for j in range(len(dense1Weights[i])):
                # uncomment the desired aggregation method

                # AVG
                #minVector.append(np.mean(np.array([dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]])))

                # MIN
                #minVector.append(min(dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]))

                # MAX
                #minVector.append(max(dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]))

                # SUM
                minVector.append(
                    dense1Weights[i][j] + dense1Weights[i+1][j] + dense1Weights[i+2][j])

        snpsDict[numSNPs] = [sbpSNPs[numSNPs], np.array(
            neuralFeatureVector), np.array(minVector)]
        numSNPs += 1

    print("========================================================================")

    numInputs = list(np.arange(int(INPUT_SIZE/3)))

    print("Initializing interaction dictionary...")
    potentialInteractions = list(combinations(numInputs, interactionWay))

    print("There are " + str(len(potentialInteractions)) +
          " potential interactions.")

    interactions = {}

    '''print("Creating empty lists...")
    for group in tqdm(potentialInteractions):
        interactions[group] = []
    print("Done.")'''


    print("Computing interaction strengths...")
    # interaction detection only for 2 way interaction (cosine similarity)
    if interactionWay == 2:

        minVecComputation = {}

        for group in potentialInteractions:
            minVecComputation[group] = []

        #------------------------#

        for intGroup in tqdm(potentialInteractions):
            minWeightsSum = 0
            for i in range(len(snpsDict[intGroup[0]][2])):
                minWeightsSum += min(snpsDict[intGroup[0]]
                                     [2][i], snpsDict[intGroup[1]][2][i])

           
            minVecComputation[intGroup] = minWeightsSum
            

            # cosine simil of neural feature vecs
            '''interactions[intGroup] = pairwise.cosine_similarity(
                (snpsDict[intGroup[0]][1]).reshape(1, -1), snpsDict[intGroup[1]][1].reshape(1, -1))'''
            # considering also min
            interactions[intGroup] = minVecComputation[intGroup]*pairwise.cosine_similarity(
                (snpsDict[intGroup[0]][1]).reshape(1, -1), snpsDict[intGroup[1]][1].reshape(1, -1))

    # interaction detection only > 2 way interaction (generalized cosine similarity for n vectors)
    if interactionWay > 2:

        minVecComputation = {}

        

        
        for intGroup in tqdm(potentialInteractions):
            
            minWeightsSum = 0
            for i in range(len(snpsDict[intGroup[0]][2])):
                minWeightsSum += min(snpsDict[intGroup[0]]
                                     [2][i], snpsDict[intGroup[1]][2][i])

            
            minVecComputation[intGroup] = minWeightsSum
            

            #####################################################################################################
            ######################################################################################################
            ########## cosine simil of neural feature vecs. Generalized version of cosine similarity ##########
            #####################################################################################################
            ######################################################################################################

            neuralFeatureVectors = []

            for i in range(interactionWay):
                neuralFeatureVectors.append(snpsDict[intGroup[i]][1])

            
            unitVectors = [(vec/np.linalg.norm(vec)).reshape(1, -1)
                           for vec in neuralFeatureVectors]

            

            concatMatrix = np.concatenate(unitVectors)

            # compute squared Frobenius norm of concatMatrix

            frobeniusNorm2 = (np.linalg.norm(concatMatrix, ord='fro'))**2

            # perform SVD on concat matrix
            u, s, vh = np.linalg.svd(concatMatrix)

            # use the generalized cosine similarity as interaction strength
            firstSigma = s[0]
            # -1 for normalization in [0,1]
            similarity = (firstSigma**2 - 1)/(frobeniusNorm2 - 1)

            
            interactions[intGroup] = minVecComputation[intGroup]*similarity
            

    # sort pairs according to interaction strenght
    if interactionWay > 3:  # unfeasible on RAM. Save file and work on disk
        intFileUnsorted = open("unsortedInteractions.txt", "w+")
        for k, v in interactions.items():
            intFileUnsorted.write(str(v) + " " + str(k) + "\n")
        intFileUnsorted.close()
        sys.exit()

    intRank = sorted(((v, k) for k, v in interactions.items()), reverse=True)

    #------------------------#
    end = time.time()
    interactingSNPs = []

    if trueData:
        if interactionWay == 2:
            for pair in intRank:
                interactingSNPs.append((pair[0], (snpsDict[pair[1][0]][0], snpsDict[pair[1][1]][0]), (
                    snpToGene[snpsDict[pair[1][0]][0]], snpToGene[snpsDict[pair[1][1]][0]])))  # add genes

        if interactionWay == 3:  # render more general for > 2
            for group in intRank:
                interactingSNPs.append((group[0], (snpsDict[group[1][0]][0], snpsDict[group[1][1]][0], snpsDict[group[1][2]][0]), (
                    snpToGene[snpsDict[group[1][0]][0]], snpToGene[snpsDict[group[1][1]][0]], snpToGene[snpsDict[group[1][2]][0]])))  # add genes

        if interactionWay == 4:  # render more general for > 2
            for group in intRank:
                interactingSNPs.append((group[0], (snpsDict[group[1][0]][0], snpsDict[group[1][1]][0], snpsDict[group[1][2]][0], snpsDict[group[1][3]][0]), (
                    snpToGene[snpsDict[group[1][0]][0]], snpToGene[snpsDict[group[1][1]][0]], snpToGene[snpsDict[group[1][2]][0]], snpToGene[snpsDict[group[1][3]][0]])))
    else:
        if interactionWay == 2:
            for pair in intRank:
                interactingSNPs.append(
                    (pair[0], (snpsDict[pair[1][0]][0], snpsDict[pair[1][1]][0])))  # add genes

        if interactionWay == 3:  # render more general for > 2
            for group in intRank:
                interactingSNPs.append(
                    (group[0], (snpsDict[group[1][0]][0], snpsDict[group[1][1]][0], snpsDict[group[1][2]][0])))  # add genes

        if interactionWay == 4:  # render more general for > 2
            for group in intRank:
                interactingSNPs.append((group[0], (snpsDict[group[1][0]][0], snpsDict[group[1][1]]
                                                   [0], snpsDict[group[1][2]][0], snpsDict[group[1][3]][0])))  # add genes

    # analyze
    print("Time elapsed for EpiCID: " + str(end-start))
    print("========================================================================")
    if save == True:
        counter = 0
        fileName = None
        if trueData:
            fileName = "../data/results/epistaticInteraction/EpiCID/BP/" + \
                str("proper_format") + \
                "/epistaticInteractions_EpiCID_" + DATASET_NAME + "_run{}.txt"
        else:
            fileName = "../data/results/epistaticInteraction/EpiCID/GAMETES/" + \
                str(DATASET_TYPE) + \
                "/epistaticInteractions_EpiCID_" + DATASET_NAME + "_run{}.txt"
        while os.path.isfile(fileName.format(counter)):
            counter += 1
        fileName = fileName.format(counter)
        intFile = open(fileName, "w+")

        for pair in interactingSNPs:

            # proper format
            if trueData:
                snp1 = str(pair[1][0])
                snp2 = str(pair[1][1])
                gene1 = str(pair[2][0])
                gene2 = str(pair[2][1])
                score = str(pair[0])[2:-2]

                intFile.write(
                    snp1 + "\t" + snp2 + "\t" + gene1 + "\t" + gene2 + "\t" + score + "\n")
            else:
                snp1 = str(pair[1][0])
                snp2 = str(pair[1][1])
                score = str(pair[0])[2:-2]
                intFile.write(snp1 + "\t" + snp2 + "\t" + score + "\n")

        print("File with interactions saved.")

    else:
        # for pair in interactingSNPs[0:100]:
        #     if "None" not in pair[2]:
        #         print(str(pair[1]) + " " + str(pair[2]) +
        #               " " + str(pair[0]) + "\n")
        pass

    print("Everything done.")

    sys.exit()

