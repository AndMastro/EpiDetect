from __future__ import absolute_import, division, print_function

import pathlib
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import random
import sys

import numpy as np
import tensorflow as tf

from tensorflow.keras import layers, Sequential, callbacks
import time

PICKLE_DIR_PATH = "../data/largeFiles/datasets/pickles/BP/" #data not present in repo since not public
DATASET_NAME = "SNPS_SBP_AVG"
PICKLE_PATH = PICKLE_DIR_PATH + DATASET_NAME + ".p"

dataset = {}
EPOCHS = 40
LEARNING_RATE = 1e-4
BATCH_SIZE = 16
INPUT_SIZE = 849  # 264*3 SBP 804 SBP, 1026 DBP, 849 PP
numLayers = "2"
save = False


def build_model():

    model = Sequential([
        layers.Dense(200, activation=tf.nn.relu, input_shape=[
                     INPUT_SIZE]),  
        layers.Dense(50, activation=tf.nn.relu),  
        #layers.Dense(1, activation=tf.nn.relu),
        layers.Dropout(0.3, noise_shape=None, seed=None),
        layers.Dense(1)
    ])

    optimizer = tf.keras.optimizers.Adam(learning_rate=LEARNING_RATE)

    model.compile(loss='mean_squared_error',  
                  optimizer=optimizer,
                  metrics=['mean_absolute_error', 'mean_squared_error'])
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

    for ind in dataset:
        xGeno = dataset[ind][0]
        if isinstance(xGeno, list):
            
            xList = xGeno
            
            trainX.append(xList)
            trainY.append(dataset[ind][2])

    print(len(trainX))

    featureAndLabels = list(zip(trainX, trainY))
    random.shuffle(featureAndLabels)
    trainX, trainY = zip(*featureAndLabels)

    X_train = np.array(trainX)
    Y_train = np.array(trainY)

    
    es = callbacks.EarlyStopping(
        monitor='val_loss', mode='min', verbose=1, patience=5, restore_best_weights=True)

    start = time.time()
    history = model.fit(
        X_train, Y_train,
        epochs=EPOCHS, validation_split=0.3, shuffle=True, verbose=1, batch_size=BATCH_SIZE)

    end = time.time()
    #plot_history(history)

    weightsLayer0 = model.layers[0].get_weights()
    weightsLayer1 = model.layers[1].get_weights()
    weightsLayer2 = model.layers[2].get_weights()
    weightsLayer3 = model.layers[3].get_weights()


    

    if save:
        np.save('../data/weights/BP/layer_0_weights_' + DATASET_NAME + '_numLayers' + numLayers,
                weightsLayer0, allow_pickle=True, fix_imports=True)
        np.save('../data/weights/BP/layer_1_weights_' + DATASET_NAME + '_numLayers' + numLayers,
                weightsLayer1, allow_pickle=True, fix_imports=True)
        np.save('../data/weights/BP/layer_2_weights_' + DATASET_NAME + '_numLayers' + numLayers,
                weightsLayer2, allow_pickle=True, fix_imports=True)
        np.save('../data/weights/BP/layer_3_weights_' + DATASET_NAME + '_numLayers' + numLayers,
                weightsLayer3, allow_pickle=True, fix_imports=True)

        print("Weights saved")

    print("Done.")
    print("Time elapsed for neural network training: ", end - start)
