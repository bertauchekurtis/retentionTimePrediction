# Predicting Retention Times of PTM Peptides using Machine Learning Techniques
### A collection of Machine Learning models for predicting the retention times of phosphorylated peptides
### Abstract
Accurately identifying protein post-translation modifications (PTMs) is important in studying cell biology
and diseases. Current experiment techniques in studying proteins and identifying peptides generate
massive amounts of data that can be extremely difficult to interpret and understand. However, it is
possible to record additional information about the protein, besides its composition, during the experiment
that can make it easier to accurately and correctly identify the protein. Specifically, the retention time of a
peptide can be recorded. Combining this measurement and comparing it to a theoretical predicted value
for the retention time can increase confidence in accurate identification.(Moruz et al., 2010) Therefore, it
is vital that there are accurate methods for predicting the retention time. Here, we explore the viability of
various types of models for predicting retention times for peptides with an emphasis on peptides with
phosphorylation modifications.
### Data
The dataset used here was created from a synthetic proteomic and phosphoproteomic dataset containing >100,000 peptides [1].
### Models
- Multiple Linear Regression
- Stepwise Linear Regression
- Linear Regression with Lasso Penalty
- Linear Regression with Ridge Penalty
- Linear Regression with Lasso and Ridge Penalties
- Random Forest
- Extreme Gradient Boostin (XG Boost)
- Support Vector Regression
#### References
[1] Marx, H., Lemeer, S., Schliep, J.E., Matheron, L., Mohammed, S., Cox, J., Mann, M., Heck, A.J. and Kuster, B., 2013. A large synthetic peptide and phosphopeptide reference library for mass spectrometry–based proteomics. Nature biotechnology, 31(6), pp.557-564.

[2] James, G., Witten, D., Hastie, T. and Tibshirani, R., 2013. An introduction to statistical learning (Vol. 112, p. 18). New York: springer.

[3] Dalpiaz, D., 2017. R for Statistical Learning.
