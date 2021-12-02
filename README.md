# retentionTimePrediction
### A collection of Machine Learning models for predicting the retetion times of phosphorylated peptides
### Data

### Linear-Based Models
- Linear Regression
    - A simple linear regression. (linearRegression.r)
- Stepwise Regression
   - A simple stepwise regression. (linearRegression.r)
- Ridge Regression
  - Penalized Regression method. (linearRegression.r)
- Lasso Regression
  - Penalized Regression method. (linearRegression.r)
- Elastic Net Regression
  - Penalized Regression method. (linearRegression.r)
### Tree-Based Models
- Random Forest
  - Tuning was performed to find best combination of number of variables and number of trees. Although a Random Forest model can be accurately assessed using Out of Bag Mean Square Error while training using the entire data set, the model was only trained with the training set to ensure comparability to other methods. (randomForest.r) The results of this tuning was visualized in a 3D plot, created in rfResults.r
- Extreme Gradient Boosting
  - Tuning was performed to find best combination of hyperparameters in xgBoost.r After the first ~2600 models, the search space was narrowed by fixing subsample, eta, and depth, all of which remained constant for the 20 best models that had been produced up to that point. Tuning continued in xgBoostHypertuningReduced.r
### Other
- Support Vector Regresion
  - Support Vector Regression was performed and the cost parameter was tuned to find the best model in the search space. This was implemented in parallel in R to help reduce computation time.
#### References
