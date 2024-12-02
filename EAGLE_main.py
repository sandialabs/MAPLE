#!/usr/bin/env python

np.random.seed(123)  # for reproducibility
import argparse
from . import loaddeepdata2
import pandas as pd
import sklearn
import scipy
from sklearn.model_selection import train_test_split
import numpy as np
import tensorflow as tf
from tensorflow import keras
from keras.layers import MultiHeadAttention, GlobalAveragePooling1D
from keras.layers import Dense, Dropout, Activation
from keras.layers import Conv1D, MaxPooling1D, Input, BatchNormalization
from keras import Model
from scikeras.wrappers import KerasClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import make_scorer, roc_auc_score


def create_CNN_with_attention(nfeats1=30, nfeats2=50, nfeats3=50, neurons1=200, neurons2=100, activation='sigmoid', reg_rate=.01, n_outputs=2, filtsize = (3,3), poolsize = 2, optimizer='Adam', dropout=0.2, init_mode='glorot_normal'):
    
    x=Input(shape=input_shape)
    p=Conv1D(nfeats1, filtsize[0], padding='same', data_format='channels_last',input_shape=input_shape,
                    kernel_initializer=init_mode)(x)
    p=BatchNormalization()(p)    
    p=Activation(activation)(p)

    ## -repeat conv and pool layers -##
    p=Conv1D(nfeats2, filtsize[1], padding='same', data_format='channels_last', kernel_initializer=init_mode)(p)
    p=BatchNormalization()(p)    
    p=Activation(activation)(p)
    p=MaxPooling1D(pool_size=poolsize, data_format='channels_last')(p)
    ##
    
    p=Conv1D(nfeats3, filtsize[1], padding='same', data_format='channels_last', kernel_initializer=init_mode,
             name='lastcnn')(p)
    p=BatchNormalization()(p)    
    p=Activation(activation)(p)
    p=MaxPooling1D(pool_size=poolsize, data_format='channels_last')(p)
    p=Activation('sigmoid')(p)

    attention_layer = MultiHeadAttention(num_heads=10, key_dim=100, name='attends')(query = p, value = p)
    d1=Dense(neurons1, trainable=True, activation=activation)(attention_layer)
    d1=Dropout(dropout)(d1)
    d2=Dense(neurons2, trainable=True, activation=activation)(d1)
    d2=Dropout(dropout)(d2)
    d2=GlobalAveragePooling1D()(d2)
    outputs=Dense(n_outputs, trainable=True, activation='softmax')(d2)
    
    model=Model(x,outputs)

    optimizer = keras.optimizers.Adam(learning_rate=1e-05)
    #compile model using metric to measure model performance
    model.compile(
        optimizer=optimizer,
        loss="categorical_crossentropy",  #"sparse_categorical_crossentropy"
        metrics=['accuracy'] #keras.metrics.AUC(multi_class = 'ovo')
        )
    return model
##

### -- HPO --- ###
def hpo_params(n_outputs=2):
    # define a grid of the hyperparameter search space ## some go into inputs above, all go into dict below
    learnRate = [1e-3, 1e-4, 1e-5]
    dropout = [0.2, 0.3, 0.4, 0.5]
    batch_size = [100,500,6000]
    epochs = [70]
    activation = ['softmax', 'softplus', 'softsign', 'relu', 'sigmoid', 'hard_sigmoid', 'linear']
    optimizer = ['SGD', 'Adagrad', 'Adam', 'Adamax', 'Nadam'] #not a model input
    momentum = [0.0, 0.25, 0.5, 0.75, 1.0]
    #init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
    init_mode = ['glorot_normal']
    reg_rate=[.1, .01, 0.001]
    neurons1 = [100,150,200,250,300]
    neurons2 = [25,50,75,100]
    nfeats1 = [15,20,25,30]
    nfeats2 = [30,40,50,60] #[10,15,20,25]
    nfeats3 = [50,75,100]

# create a dictionary from the hyperparameter grid 
    grid = dict(
        model__dropout=dropout,
        batch_size=batch_size,
        epochs=epochs,
        optimizer=optimizer,
        optimizer__learning_rate=learnRate,
        optimizer__momentum=momentum,
        model__init_mode=init_mode,
        model__activation=activation,
        model__reg_rate=reg_rate,
        model__neurons1 = neurons1,
        model__neurons2 = neurons2,
        model__nfeats1= nfeats1,
        model__nfeats2=nfeats2,
        model__nfeats3=nfeats3
    )
    # optimizer = keras.optimizers.Adam(learning_rate=0.001)
    # custom metrics
    if n_outputs == 2:
        auc_score = make_scorer(roc_auc_score, needs_proba=True)
                             
    elif n_outputs > 2:
        auc_score = make_scorer(roc_auc_score,
                            #needs_threshold=True,
                            needs_proba=True,
                            multi_class = 'ovr', average='micro') 

    return auc_score, grid
#----HPO end-----#


##Default values based on hyperparameter optimization
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("species1", default="Ncrassa", type=str)
    parser.add_argument("species2", default="Fgram", type=str)
    parser.add_argument("n_outputs", default=2, type=int)
    parser.add_argument("bin_num", default=50, type=int)
    parser.add_argument("--species1cov", default=None)
    parser.add_argument("--species1rna", default=None)
    parser.add_argument("--species1mods", default=None, type=list)
    parser.add_argument("--species2cov", default=None)
    parser.add_argument("--species2rna", default=None)
    parser.add_argument("--species2mods", default=None, type=list)

    parser.add_argument("--epochs", default=100, type=int)
    parser.add_argument("--batch_size", default=100, type=int)
    parser.add_argument("--trained_model_result", type=str)
    parser.add_argument("--runHPO", default=False, type=bool)

    args = parser.parse_args()

    # X_train1, X_test1, y_train1, y_test1, X_train2, X_test2, y_train2, y_test2, savegene_train1,savegene_test1, savegene_train2,savegene_test2 = (loaddeepdata2_wGenes4GO.loaddeepdata('LmacL','Ncrassa',2))
    X_train1, X_test1, y_train1, y_test1, X_train2, X_test2, y_train2, y_test2 = loaddeepdata2.loaddeepdata(args.species1, args.species2,  
                                                                                                            args.n_outputs, args.bin_num, 
                                                                                                            args.species1cov, args.species1rna, args.species1mods, 
                                                                                                            args.species2cov, args.species2rna, args.species2mods)

    input_shape = (X_train1.shape[1], X_train1.shape[2]) # model_attention = create_CNN_with_attention()
    n_outputs = args.n_outputs
    model = create_CNN_with_attention(n_outputs=n_outputs)
    print(model.summary())

    if args.runHPO is True:
        auc_score, grid = hpo_params(n_outputs=n_outputs)
        hpo_model = KerasClassifier(model=model, verbose=0)

        # initialize a random search with a 3-fold cross-validation and then
        # start the hyperparameter search process
        print("[INFO] performing random search...")
        searcher = RandomizedSearchCV(estimator=model,param_distributions=grid,error_score="raise", cv=10,
                                    scoring=auc_score, return_train_score=True)

        searchResults = searcher.fit(X_train1, y_train1, validation_data=(X_test1, y_test1))

        # summarize grid search information
        bestScore = searchResults.best_score_
        bestParams = searchResults.best_params_
        print("[INFO] best score is {:.2f} using {}".format(bestScore, bestParams))

        print("[INFO] evaluating the best model...")
        model = searchResults.best_estimator_
        accuracy = model.score(X_test2, y_test2)
        print("accuracy: {:.2f}%".format(accuracy * 100))

    else:
        model = model.fit(X_train1, y_train1, validation_data=(X_test1, y_test1), 
                            epochs=args.epochs, batch_size=args.batch_size)

    if args.trained_model_result is not None:
        model.save(f'{args.trained_model_result}.keras')
        
    #Metrics considered in study
    # accuracy = model.evaluate(X_test2, y_test2)
    # #print("accuracy: {:.2f}%".format(accuracy * 100))
    # ypred =  model.predict(X_test2)
    # sklearn.metrics.roc_auc_score(y_test2, ypred)#, average='micro', multi_class='ovr')
