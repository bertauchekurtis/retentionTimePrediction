# Predicting Retention Times of PTM Peptides using Machine Learning Techniques
### A collection of models for predicting the retention times of peptides, with a focus on post-translational phosphorylation modifications
### Abstract
Accurately identifying protein post-translation modifications (PTMs) is important in studying cell biology
and diseases. Current experiment techniques in studying proteins and identifying peptides generate
massive amounts of data that can be extremely difficult to interpret and understand. However, it is
possible to record additional information about the protein, besides its composition, during the experiment
that can make it easier to accurately and correctly identify the protein. Specifically, the retention time of a
peptide can be recorded. Combining this measurement and comparing it to a theoretical predicted value
for the retention time can increase confidence in accurate identification. (Moruz et al., 2010) Therefore, it
is vital that there are accurate methods for predicting the retention time. Here, we explore the viability of
various types of models for predicting retention times for peptides with an emphasis on peptides with
phosphorylation modifications.

### Models
- Multiple Linear Regression
- Stepwise Linear Regression
- Linear Regression with Lasso Penalty
- Linear Regression with Ridge Penalty
- Linear Regression with Lasso and Ridge Penalties
- Random Forest
- Extreme Gradient Boostin (XG Boost)
- Support Vector Regression
- Bayesian Linear Regression

### Files
Each model is contained within its own script. Inside of each script is both the training and tuning procedure for each model as well as the evalutation for each model. In addtion to the model files, there are also files for other miscellaneous tasks, such as splitting the data and creating other file formats.

### References

Bouwmeester, R., Gabriels, R., Hulstaert, N., et al. (2021). Deeplc can predict retention times for peptides188
that carry as-yet unseen modifications. Nat Methods, 18:1363–1369.189

Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., Cho, H., Chen, K., Mitchell, R., Cano, I., Zhou,190
T., Li, M., Xie, J., Lin, M., Geng, Y., and Li, Y. (2021). xgboost: Extreme Gradient Boosting. R191
package version 1.4.1.1.192

Creasy, D. M. and Cottrell, J. S. (2004). Unimod: Protein modifications for mass spectrometry. PRO-193
TEOMICS, 4(6):1534–1536.194

Dalpiaz, D. (2020). R for statistical learning. https://daviddalpiaz.github.io/r4sl/t.195

Gessulat, S., Schmidt, T., Zolg, D. P., et al. (2019). Prosit: proteome-wide prediction of peptide tandem196
mass spectra by deep learning. Nat Methods, 16:509–518.197

Guan, S., Moran, M. F., and Ma, B. (2019). Prediction of lc-ms/ms properties of peptides from sequence198
by deep learning. Molecular cellular proteomics : MCP, 18(10):2099–2107.199

James, G., Witten, D., Hastie, T., and Tibshirani, R. (2013). An Introduction to Statistical Learning.200
Springer, New York.201

Lumley, T. and based on Fortran code by Alan Miller (2020). leaps: Regression Subset Selection. R202
package version 3.1.203

Marx, H., Lemeer, S., Schliep, J., et al. (2013). A large synthetic peptide and phosphopeptide reference204
library for mass spectrometry–based proteomics. Nat Biotechnol, 31:557–564.205

Meek, J. L. (1980). Prediction of peptide retention times in high-pressure liquid chromatography on the206
basis of amino acid composition. Proceedings of the National Academy of Sciences, 77(3):1632–1636.207

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., and Leisch, F. (2021). e1071: Misc Functions208
of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien. R package209
version 1.7-9.210

Moruz, L. and K  ̈all, L. (2017). Peptide retention time prediction. Mass Spectrometry Reviews, 36(5):615–211
623.212

Moruz, L., Tomazela, D., and K  ̈all, L. (2010). Training, selection, and robust calibration of retention time213
models for targeted proteomics. Journal of Proteome Research, 9(10):5209–5216. PMID: 20735070.214
5/6

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer,215
P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., and216
Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. Journal of Machine Learning217
Research, 12:2825–2830.218

R Core Team (2021). R: A Language and Environment for Statistical Computing. R Foundation for219
Statistical Computing, Vienna, Austria.220

Sch  ̈onthal, A. H. (2001). Role of serine/threonine protein phosphatase 2a in cancer. Cancer Letters,221
170(1):1–13.222

Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011). Regularization paths for cox’s proportional223
hazards model via coordinate descent. Journal of Statistical Software, 39(5):1–13.224

Wen, B., Li, K., Zhang, Y., et al. (2020). Cancer neoantigen prioritization through sensitive and reliable225
proteogenomics analysis. Nat Commun, 11:1759.226

Wright, M. N. and Ziegler, A. (2017). ranger: A fast implementation of random forests for high227
dimensional data in C++ and R. Journal of Statistical Software, 77(1):1–17.228

Yan, J., Packer, N., Gooley, A., and Williams, K. (1998). Protein phosphorylation: technologies for the229
identification of phosphoamino acids. Journal of Chromatography A, 808(1-2):23–41.
