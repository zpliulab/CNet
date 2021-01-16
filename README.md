0. network_match.R ---- Network structure of DEGs in RegNetwork
1. Data_net.R ---- Intersection of genes in the Data and genes in RegNetwork
2. DE_net.R ---- Select differentially expressed genes (DEGs) and consider the network structure in RegNetwork
3. data_split2.R ----  Randomly divide Dataset into 20 new Datasets by 1:1 of the positive and negative samples
4. cut.R ---- Diameters and cut-nodes of component of DEGs in RegNetwork
5. adj.R ---- Adjacency matrix and its eigenvalues
6. CNet_RLR_126e50102.m ---- Main function
7. coef2feature_data ---- Select feature genes based on threshold and Extract feature expression values in test dataset and independent dataset
8. class_net.R ---- Independent data set verification to obtain ROC curve and classification accuracy

costFunction12.m ---- Objective function

cv.m ---- Cross validation to select optimal parameters

Laplacian_Matrix.m ---- Laplacian matrix according to the adjacency matrix

LogitisLap.m ---- LogitisLap function for cv

SGNLR.m ---- SGNLR function for LogitisLap

my_error.m ---- error function 

getLambMax.m ---- getLambMax function for cv

Predict.m ---- Predict function on test dataset

plot_roc.m ---- Roc curve function on test dataset

RunRcode.m ---- RunRcode for main function
