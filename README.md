# Thesis
Code for *Machine learning stractification for oncology patient survival*.

This code generates/extracts the data and plots for the thesis. 
All code is written in R.

To run, all files are required.
In each chapter, each plot has an associated script:
- runCh1Fig3.R
- runCh3Fig1.R
- runCh3Fig2.R
- runCh3Fig3.R
- runCh4SyntheticExp1.R
- runCh4SyntheticExp2.R
- runCh4SyntheticExp3.R
- runCh4RealTothillEtAl.R
- runCh4RealYuanEtAl.R
- runCh5SyntheticExp1IR.R
- runCh5SyntheticExp2R.R
- runCh5SyntheticExp3R.R
- runCh5RealTothillEtAl.R

These scripts may be found in the 'scripts' folder. 
Any necessary functions will be sourced from the 'toSource' folder, which is required to be in the same folder that contains the working directory.

Upon completion, each script will save results into a folder named "Runs" in the working directory, which is required to already exist.

Both Gaussian process regression and the Gaussian process for survival data models, GPS1, GPS2 and GPS3, may be run using the ApplyGP.R function, as applied in the scripts above. The feature selection methods IARD, RSFS and GPSurvBIC also use this function.

For the real data experiments, data and/or code were used as supplied by the original authors where possible. 
For the Yuan et al. (2014) data, code was obtained from Synapse, synapse id:syn1720423. Data is extracted from the same source using the Synapse R client. This requires a login.
The Tothill et al. (2008) data was extracted using the R package curatedOvarianData, from the data set GSE9891.

Required libraries are:
- akima
- caret
- curatedOvarianData
- devtools
- fields
- foreach
- gbm
- glmnet
- impute
- ipred
- MASS
- Matrix
- mice
- nlme
- NORMT3
- pdist
- randomForestSRC
- rgl
- rms
- survAUC
- survcomp
- survival
- synapseClient
- UpSetR
- zoo
