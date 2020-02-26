import numpy as np
import tempfile
import os
from itertools import combinations
import itertools as IT
from snpReader import snpRead
from snpReader import geneRead
from sklearn.metrics import pairwise
import statistics
import sys
import math
from operator import add
import matplotlib.pyplot as plt
from tqdm import tqdm

INPUT_SIZE = 256*3  # 804 SBP, 1026 DBP, 849 PP
DATASET_NAME = "threshold_risk_model_MAF01_eta01_theta0.5_EDM-1_08"
if __name__ == "__main__":

    trueData = False  # variable stating if using true or generated data

    # if u wanna save on file
    save = False

    if sys.argv[1] == "save":
        save = True

    # specify degree of interaction

    interactionWay = int(sys.argv[2])

    numLayers = int(sys.argv[3])

    layerWeights0 = None
    layerWeights1 = None
    layerWeights2 = None
    layerWeights3 = None
    layerWeights4 = None
    layerWeights5 = None
    sbpSNPs = None

    layerWeights0 = np.load(
        '..\\data\\weights\\layer_0_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)
    layerWeights1 = np.load(
        '..\\data\\weights\\layer_1_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)
    layerWeights2 = np.load(
        '..\\data\\weights\\layer_2_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)
    layerWeights3 = np.load(
        '..\\data\\weights\\layer_3_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)

    if numLayers == 4:
        layerWeights4 = np.load(
            '..\\data\\weights\\layer_4_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)
        layerWeights5 = np.load(
            '..\\data\\weights\\layer_5_weights_' + DATASET_NAME + '_numLayers' + str(numLayers) + '.npy', allow_pickle=True)

    print("Loading SNPs database...")

    sbpSNPs = snpRead(
        "..\\data\\snpList256.bim")

    dense1 = None
    dense2 = None
    dense3 = None
    dense4 = None
    output = None

    if numLayers == 4:
        dense1 = layerWeights0[0]
        dense2 = layerWeights1[0]
        dense3 = layerWeights2[0]
        dense4 = layerWeights3[0]
        #layerWeright4 is dropout (empty)
        output = layerWeights5[0]
    else:
        dense1 = layerWeights0[0]
        dense2 = layerWeights1[0]
        #layerWeright2 is dropout (empty)
        output = layerWeights3[0]

    print("===========================")

    print("Done!")
    print("===========================")

    # Define weights vectors NOT according to NID paper
    # NOT LIKE THAT: we need all the weights for all the input variables for any unit. So we need a matrix such thta we have a vector of 153 elems for any of the 64 hidden units
    #Matrix is not transposed
    print("Defining vectors for interaction detection. We need absolute values...")

    dense1Weights = []
    firstAppend = True

    # for any input unit we have a vectyro with the weights to all the fist layer hidden units
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

    if numLayers == 4:
        dense3Weights = []
        firstAppend = True

        for inputVec in dense3:
            for i in range(len(inputVec)):
                if firstAppend:
                    dense3Weights.append([abs(inputVec[i])])
                else:
                    dense3Weights[i].append(abs(inputVec[i]))
            firstAppend = False

        dense4Weights = []
        firstAppend = True

        for inputVec in dense4:
            for i in range(len(inputVec)):
                if firstAppend:
                    dense4Weights.append([abs(inputVec[i])])
                else:
                    dense4Weights[i].append(abs(inputVec[i]))
            firstAppend = False

    outputWeights = []
    for w in output:
        outputWeights.append((list(abs(w))))

    def flatten(l): return [item for sublist in l for item in sublist]

    outputWeights = flatten(outputWeights)

    print("Defining Aggregated Weight for any unit...")

    aggregateWeightsD1 = None

    if numLayers == 4:
        aggregateWeightsD1 = np.matmul(
            np.array(outputWeights), np.array(dense4Weights))
        aggregateWeightsD1 = np.matmul(
            aggregateWeightsD1, np.array(dense3Weights))
        aggregateWeightsD1 = np.matmul(
            aggregateWeightsD1, np.array(dense2Weights))
    else:
        aggregateWeightsD1 = np.matmul(
            np.array(outputWeights), np.array(dense2Weights))

    print("Done.")

    # defining aggregated neural feature vector for any SNP (aggregation of 3 input features as follows)
    numSNPs = 0
    snpsDict = {}
    scalingFactor11 = 0.5
    scalingFactor12 = 1
    scalingFactor22 = 2
    for i in range(0, len(dense1Weights), 3):
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
            neuralFeatureVector, neuralFeatureVectorProxy3*scalingFactor22)

        # uncomment this line to have NFV concatenation
        #neuralFeatureVector = np.concatenate((neuralFeatureVectorProxy1, neuralFeatureVectorProxy2, neuralFeatureVectorProxy3*scalingFactor22))

        # uncomment this to have max
        #neuralFeatureVector = []
        '''for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(max(neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]*scalingFactor22))
        minVector = []'''

        # for min
        '''neuralFeatureVector = []
        for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(min(neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]*scalingFactor22))
        minVector = []'''

        # for avg
        '''neuralFeatureVector = []
        for i in range (0,len(neuralFeatureVectorProxy1)):
            neuralFeatureVector.append(np.mean(np.array((neuralFeatureVectorProxy1[i],neuralFeatureVectorProxy2[i],neuralFeatureVectorProxy3[i]*scalingFactor22))))
        minVector = []'''

        minVector = []
        concat = False

        if concat:
            minVector = np.concatenate(
                (dense1Weights[i], dense1Weights[i+1], dense1Weights[i+2]*scalingFactor22))
        else:
            for j in range(len(dense1Weights[i])):
                # uncomment the desired aggregation method

                # AVG
                #minVector.append(np.mean(np.array([dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]])))

                # MIN
                #minVector.append(min(dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]))

                # MAX
                #minVector.append(max(dense1Weights[i][j], dense1Weights[i+1][j], dense1Weights[i+2][j]*scalingFactor22))

                # SUM
                minVector.append(
                    dense1Weights[i][j] + dense1Weights[i+1][j] + dense1Weights[i+2][j])

        snpsDict[numSNPs] = [sbpSNPs[numSNPs], np.array(
            neuralFeatureVector), np.array(minVector)]
        numSNPs += 1

    print("========================================================================")

    numInputs = list(np.arange(int(INPUT_SIZE/3)))

    print("Inizializing interaction dictionary...")
    potentialInteractions = list(combinations(numInputs, interactionWay))

    print("There are " + str(len(potentialInteractions)) +
          " potential interactions.")

    interactions = {}

    '''print("Creating empty lists...")
    for group in tqdm(potentialInteractions):
        interactions[group] = []
    print("Done.")'''

    # the following lines of code are  used to study the interaction strength according to the min vec only to check how to exploit it.
    #------------------------#

    print("Computing interaction strengths...")
    # interaction detection only for 2 way interaction (cosine similarity)
    if interactionWay == 2:

        minVecComputation = {}

        for group in potentialInteractions:
            minVecComputation[group] = []

        #------------------------#

        # in order to exploit FUTURE (deeper layers) chek old epistatic intreaction file from Master thesis
        for intGroup in tqdm(potentialInteractions):
            minWeightsSum = 0
            for i in range(len(snpsDict[intGroup[0]][2])):
                minWeightsSum += min(snpsDict[intGroup[0]]
                                     [2][i], snpsDict[intGroup[1]][2][i])

            # the following lines of code are  used to study the interaction strength according to the min vec only to check how to exploit it.
            #------------------------#
            minVecComputation[intGroup] = minWeightsSum
            #------------------------#

            # cosine simil of neural feature vecs
            '''interactions[intGroup] = pairwise.cosine_similarity(
                (snpsDict[intGroup[0]][1]).reshape(1, -1), snpsDict[intGroup[1]][1].reshape(1, -1))'''
            # considering also min
            interactions[intGroup] = minVecComputation[intGroup]*pairwise.cosine_similarity(
                (snpsDict[intGroup[0]][1]).reshape(1, -1), snpsDict[intGroup[1]][1].reshape(1, -1))

    # interaction detection only > 2 way interaction (generalized cosine similarity for n vectors)
    if interactionWay > 2:

        #minVecComputation = {}

        '''for group in potentialInteractions:
            minVecComputation[group] = []'''

        #------------------------#

        # in order to exploit FUTURE (deeper layers) chek old epistatic intreaction file from Master thesis
        for intGroup in tqdm(potentialInteractions):
            # minweight sum not needed. So do not compute it.
            #minWeightsSum = 0
            '''for i in range(len(snpsDict[intGroup[0]][2])):
                minWeightsSum += min(snpsDict[intGroup[0]]
                                    [2][i], snpsDict[intGroup[1]][2][i])'''

            # the following lines of code are  used to study the interaction strength according to the min vec only to check how to exploit it.
            #------------------------#
            #minVecComputation[intGroup] = minWeightsSum
            #------------------------#

            #####################################################################################################
            ######################################################################################################
            ########## cosine simil of neural feature vecs. Generalized version of cosine similarity ##########
            #####################################################################################################
            ######################################################################################################

            neuralFeatureVectors = []

            for i in range(interactionWay):
                neuralFeatureVectors.append(snpsDict[intGroup[i]][1])

            # define unit vectors as vec/norm(vec)  shape is: (columns/rows/pages/blocks/...)
            unitVectors = [(vec/np.linalg.norm(vec)).reshape(1, -1)
                           for vec in neuralFeatureVectors]

            # concat vectors into a matrix

            concatMatrix = np.concatenate(unitVectors)

            # compute squared Frobenius norm of concatMatrix

            frobeniusNorm2 = (np.linalg.norm(concatMatrix, ord='fro'))**2

            # perform SVD on concat matrix
            u, s, vh = np.linalg.svd(concatMatrix)

            # use the generalized cosine similarity as interaction strength
            firstSigma = s[0]
            # -1 for normalization in [0,1]
            similarity = (firstSigma**2 - 1)/(frobeniusNorm2 - 1)

            interactions[intGroup] = similarity
            # considering also min
            #interactions[intGroup] = minVecComputation[intGroup]*pairwise.cosine_similarity((snpsDict[intGroup[0]][1]).reshape(1,-1), snpsDict[intGroup[1]][1].reshape(1,-1))

    # sort pairs according to interaction strenght
    if interactionWay > 3:  # unfeasible on RAM. Save file and work on disk
        intFileUnsorted = open("unsortedInteractions.txt", "w+")
        for k, v in interactions.items():
            intFileUnsorted.write(str(v) + " " + str(k) + "\n")
        intFileUnsorted.close()
        sys.exit()

    intRank = sorted(((v, k) for k, v in interactions.items()), reverse=True)

    #------------------------#

    interactingSNPs = []

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

    # and now, analyze

    if save == True:
        counter = 0
        fileName = "..\\data\\results\\epistaticInteraction\\synthetic\\" + \
            str(interactionWay) + \
            "\\epistaticInteractions_thesisMethodandMinVecSum_" + DATASET_NAME + "_run{}.txt"
        while os.path.isfile(fileName.format(counter)):
            counter += 1
        fileName = fileName.format(counter)
        intFile = open(fileName, "w+")

        for pair in interactingSNPs:
            intFile.write(
                str(pair[1]) + " " + str(pair[0]) + "\n")

        print("File with interactions saved.")

    else:
        for pair in interactingSNPs[0:100]:
            if "None" not in pair[2]:
                print(str(pair[1]) + " " + str(pair[2]) +
                      " " + str(pair[0]) + "\n")

    interValues = []
    for val in interactingSNPs:
        interValues.append((val[0].tolist())[0][0])

    plt.hist(interValues)

    saveLocation = "../data/results/plots/strengthDistribution" + DATASET_NAME + ".png"
    plt.savefig(saveLocation)

    plt.show()

    print("Everything done.")

    sys.exit()

    ################################################################################
    ################################################################################
    ##################              FILE ENDS HERE            ######################
    ################################################################################
    ################################################################################
