#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import time
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import cohen_kappa_score
from sklearn.preprocessing import StandardScaler
import pickle
import random

random.seed(42)
print(random.random())


# It is needed to have the class column
train = pd.read_csv("train.tsv", sep=',', index_col = 0)

print(train['class'])

X = train.drop(columns='class')
Y = train['class']
Y = Y.astype('category')
X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.30)

print(X_train.shape, X_val.shape, Y_train.shape, Y_val.shape)

scaler = StandardScaler()

normalized_X_train = pd.DataFrame(
    scaler.fit_transform(X_train),
    columns = X_train.columns
)

model = XGBClassifier(
    objective= 'binary:logistic'
)

param_grid = {
    'max_depth': range (7, 15, 20),
    'n_estimators': range(150, 300, 500),
    'learning_rate': [0.01, 0.1, 0.2, 0.3]
}

grid_search = GridSearchCV(estimator=model,
                             param_grid=param_grid,
                             cv=StratifiedKFold(5),
                             n_jobs = 6, verbose = 0,
                             scoring = 'accuracy')

grid_search.fit(normalized_X_train, Y_train)

xgb_model = grid_search.best_estimator_
xgb_model.fit(normalized_X_train,Y_train)


normalized_X_val = pd.DataFrame(
    scaler.fit_transform(X_val),
    columns = X_val.columns
)


Y_pred = xgb_model.predict(normalized_X_val)


print("Accuracy")
print(accuracy_score(Y_val, Y_pred))
print("Kappa")
print(cohen_kappa_score(Y_val, Y_pred))


pkl_filename = "xgb_model.pkl"
with open(pkl_filename, 'wb') as file:
    pickle.dump(xgb_model, file)
