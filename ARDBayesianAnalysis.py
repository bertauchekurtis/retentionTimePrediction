# libraries
import math
import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from sklearn import linear_model
from sklearn.linear_model import Ridge

from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import median_absolute_error
from sklearn.metrics import r2_score

random.seed(37)
# importing data ==================================
df = pd.read_csv("/Users/rileybrenner/Desktop/trainingSet.csv")
df_test = pd.read_csv("/Users/rileybrenner/Desktop/testingSet.csv")
print("|==========> Data Loaded <==========|")
# =================================================

#print(df.head())
#print(df_test.head())



print("\n|==========> Creating Matrix (xgboost) <==========|")
RetTime = df["RetentionTime"]
data = df.drop(["RetentionTime", "Peptide.Sequence2"], axis = 1)

RetTimeTest = df_test["RetentionTime"]
data_test = df_test.drop(["RetentionTime", "Peptide.Sequence2"], axis = 1)
print(data_test.head())

y_vals = []
for x in RetTime:
    y_vals.append(int(x))

y_vals_test = []
for x in RetTimeTest:
    y_vals_test.append(int(x))

print(data.shape)
print(RetTime)

# #############################################################################
# Fit the Bayesian Ridge Regression and an OLS for comparison
clf = linear_model.ARDRegression()

clf.fit(data, y_vals)

y_vals_preidct = clf.predict(data_test)

print("Score:", clf.score(data_test, y_vals_test))

print(y_vals_preidct)

print("Mean Absolyte Error:",mean_absolute_error(y_vals_test, y_vals_preidct))
print("Mean Squared Error:",mean_squared_error(y_vals_test, y_vals_preidct))
print("RMSE:" , math.sqrt(mean_squared_error(y_vals_test, y_vals_preidct)))
print("Median Absolute Error:",median_absolute_error(y_vals_test, y_vals_preidct))
print("R^2:",r2_score(y_vals_test, y_vals_preidct))

print(max(y_vals)) 
print(min(y_vals))

collected_errors = y_vals_preidct - y_vals_test

plt.hist(collected_errors)

q1=np.quantile(residuals, 0.025)

q2=np.quantile(residuals, 0.975)

print("Lower Qunatile:", q1)
print("Upper Quantile:", q2)
print("Window width:", abs(q1)+abs(q2))
print("Correlation:",np.corrcoef(y_vals_preidct, y_vals_test))