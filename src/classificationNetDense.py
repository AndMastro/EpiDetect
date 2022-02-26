from __future__ import absolute_import, division, print_function

import pickle
import matplotlib.pyplot as plt
import pandas as pd
import random
import sys
import numpy as np
import tensorflow as tf
import tensorflow.keras
from tensorflow import keras
from tensorflow.keras import layers
from tqdm import tqdm


DATASET_TYPE = "marginal_effect"
PICKLE_DIR_PATH = "../data/largeFiles/datasets/pickles/GAMETES/" + DATASET_TYPE + "/"
DATASET_NAME = "epistatic_plus_ME_risk_model_MAF01_eta01_theta1_lambda1_2_EDM-1_01" #1 is MAF=0.2 and 2 is MAF=0.5 ---- it should make sense now. Prevalence should reflect #controls. Am I right?
PICKLE_PATH = PICKLE_DIR_PATH + DATASET_NAME + ".p"
###RENDER DATASET INPUTABLE. IT should be more modular###

########################################
# arguments in input: numLayers (2 or 4)

dataset = {}
EPOCHS = 30
LEARNING_RATE = 1e-4
BATCH_SIZE = 64
INPUT_SIZE = 256*3 #2000*3 256*3
numLayers = sys.argv[1]


def build_model():

    model = None

    if numLayers == '4':

        # classification
        model = keras.Sequential([
            layers.Dense(140, activation=tf.nn.relu,
                         input_shape=[INPUT_SIZE], kernel_initializer='glorot_uniform'),
            layers.Dense(100, activation=tf.nn.relu, kernel_initializer='glorot_uniform'),
            layers.Dense(60, activation=tf.nn.relu, kernel_initializer='glorot_uniform'),
            layers.Dense(20, activation=tf.nn.relu, kernel_initializer='glorot_uniform'),
            layers.Dropout(0.3, noise_shape=None, seed=None),
            # layers.Dense(2, activation='softmax') #multi-class classification
            layers.Dense(1, activation='sigmoid')  # binary classification
        ])

    elif numLayers == '2':
        model = keras.Sequential([
            layers.Dense(200, activation=tf.nn.relu,
                         input_shape=[INPUT_SIZE]),
            layers.Dense(100, activation=tf.nn.relu),
            layers.Dropout(0.3, noise_shape=None, seed=None),
            # layers.Dense(2, activation='softmax')  #multi-class classification
            layers.Dense(1, activation='sigmoid')  # binary classification
        ])

    else:
        print("Command not recognized. Aborting...")
        sys.exit()

    optimizer = tf.keras.optimizers.Adam(learning_rate=LEARNING_RATE)
    #optimizer = tf.keras.optimizers.Adam()
    #optimizer = tf.keras.optimizers.SGD()

    # classification
    model.compile(loss='binary_crossentropy',  # sparse_categorical_crossentropy binary_crossentropy
                  optimizer=optimizer,
                  metrics=['accuracy'])

    return model


def plot_history(history):
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch

    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Abs Error')
    plt.plot(hist['epoch'], hist['mean_absolute_error'],
             label='Train Error')
    plt.plot(hist['epoch'], hist['val_mean_absolute_error'],
             label='Val Error')
    # plt.ylim([0,0.125])
    plt.ylim([0, 25])
    plt.legend()

    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Square Error')
    plt.plot(hist['epoch'], hist['mean_squared_error'],
             label='Train Error')
    plt.plot(hist['epoch'], hist['val_mean_squared_error'],
             label='Val Error')
    # plt.ylim([0,0.025])
    plt.ylim([0, 625])
    plt.legend()

    plt.show()


if __name__ == "__main__":
    model = build_model()
    trainX = []
    trainY = []

    with open(PICKLE_PATH, 'rb') as handle:
        dataset = pickle.load(handle)

    print("Loading dataset...")

    for ind in tqdm(dataset):
        #xGeno = [1,0,0]*256
        xGeno = dataset[ind][0]  # one-hot encoding for the SNPs
        #xGeno += dataset[ind][762:]
        trainX.append(xGeno)
        label = dataset[ind][1]  #index is 2 if we use real SBP, PP andDBP data (1 is for pheno). Index is 1 for synthetic data.
        trainY.append(label)
    print("Dataset loaded. Num samples: ")
    print(len(trainX))

    #print("Balancing dataset...")
    # randominze data to prevent bias
    featureAndLabels = list(zip(trainX, trainY))
    random.shuffle(featureAndLabels)
    trainX, trainY = zip(*featureAndLabels)

    '''X_train, X_test, y_train, y_test = train_test_split(
        trainX, trainY, test_size=0.3, random_state=42)'''

    # uncomment to balance dataset
    '''trainX_balanced = []
    trainY_balanced = 
    minSamples = min(numSamples)
    print("Each class has samples: ")
    print(numSamples)'''

    '''taken = [0]*numClasses
    for i in tqdm(range(len(trainX))):
        if taken[trainY[i]] >= minSamples:
            continue
        trainX_balanced.append(trainX[i])
        trainY_balanced.append(trainY[i])
        taken[trainY[i]] += 1

    print("Done. Num sample per class: " + str(minSamples))'''

    # cast from list to numpy arrays. Balanced version
    '''X_train = np.array(trainX_balanced)
    Y_train = np.array(trainY_balanced)'''

    # unbalanced version
    X_train = np.array(trainX)
    y_train = np.array(trainY)
    #X_test = np.array(X_test)
    #y_test = np.array(y_test)

    #X_train = np.reshape(X_train, (X_train.shape[0], -1))

    '''print("DEBUG PRINT:")
    print(X_train[0])
    print(Y_train[0])

    sys.exit()'''

    # define early stopping callback
    es = keras.callbacks.EarlyStopping(
        monitor='val_accuracy', mode='max', verbose=1, patience=5, restore_best_weights=True) #val-accuracy or val_loss? USe acc up to know, so be coherent. Frot he future, loss

    # set class weights in order to work with unbalanced data. We keep thus statistical power
    '''print("Computing class weigth to work with unbalanced data...")
    class_weights = class_weight.compute_class_weight('balanced',
                                                      np.unique(Y_train),
                                                      Y_train)
    print("Class weights:")
    print(class_weights)'''

    history = model.fit(
        X_train, y_train,
        epochs=EPOCHS, validation_split = 0.3, shuffle=True, verbose=1, batch_size=BATCH_SIZE, callbacks=[es])

    # plot_history(history)

    weightsLayer0 = None
    weightsLayer1 = None
    weightsLayer2 = None
    weightsLayer3 = None
    weightsLayer4 = None
    weightsLayer5 = None

    if numLayers == '2':
        weightsLayer0 = model.layers[0].get_weights()
        weightsLayer1 = model.layers[1].get_weights()
        weightsLayer2 = model.layers[2].get_weights()
        weightsLayer3 = model.layers[3].get_weights()

    else:
        weightsLayer0 = model.layers[0].get_weights()
        weightsLayer1 = model.layers[1].get_weights()
        weightsLayer2 = model.layers[2].get_weights()
        weightsLayer3 = model.layers[3].get_weights()
        weightsLayer4 = model.layers[4].get_weights()
        weightsLayer5 = model.layers[5].get_weights()

    save = True

    if save:
        if numLayers == '4':
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_0_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer0, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_1_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer1, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_2_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer2, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_3_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer3, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_4_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer4, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_5_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer5, allow_pickle=True, fix_imports=True)
        else:
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_0_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer0, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_1_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer1, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_2_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer2, allow_pickle=True, fix_imports=True)
            np.save('../data/weights/GAMETES/' + DATASET_TYPE + '/layer_3_weights_'+ DATASET_NAME + '_numLayers' + numLayers,
                    weightsLayer3, allow_pickle=True, fix_imports=True)

        print("Weights saved")

    print("Done.")
