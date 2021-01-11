rm(list=ls())

library(pROC)
library(MASS)
library(glmnet)
# install.packages("glmnet")
# ------------------------------------------
# Please first set the path
working_dir = "C:/Users/LiLingyu/Desktop/NSLR-master"
setwd(working_dir)

source('Fun_Auxiliary.R')
source('Fun_NSLR.R')

# ------------------------------------------
# Generate a simple simulation data
snrlam0 = 3
f.num0  = 100
sim.data = GET.SIM.DATA2(smaple.num = 700, feature.num=f.num0, random.seed = 10, snrlam=snrlam0)
# adj = get_sim_prior_Net(f.num0, 40, 0.3, 0.05)
adj = get_sim_prior_Net(f.num0, 40, 0.05, 0.02)
# write.csv(adj, file = "adj_example.csv", row.names = F)

Train.id = 1:300
Valid.id = 301:500
w.true = sim.data$w

# Training data
X1 = sim.data$X[Train.id,]; y1 =sim.data$y[Train.id]; 
train_data <- as.matrix(cbind(y1, X1))
colnames(train_data) <- c("Label", c(1:100))
dim(train_data)
# View(train_data[,1:10])
# write.csv(train_data, file = "Data_train\\1.txt")


# Testing data
X2 = sim.data$X[Valid.id,]; y2 =sim.data$y[Valid.id]; 
test_data <- cbind(y2,X2)
colnames(test_data) <- c("Label", c(1:100))
dim(test_data)
# View(test_data[,1])
# write.csv(test_data, file = "Data_test\\1.txt")


# Non Normalized Laplacian Matrix"
PM = get.penalityMatrix(adj,X1, y1)

# ----------------------------------------
# ----------------------------------------
# Two regularization parameters
lambda = 0.65
alpha  = 0.3

# ----------------------------------------
# A typical logistic regression model and four regularized logistic regession models
# ----------------------------------------

# Typical Logistic Regression Model
out1 =    SGNLR(X1, y1, PM$M.c,        lambda=0, alpha=0, niter=20)

# L1 (Lasso) Regularized Logistic Regression Model
out2 =    SGNLR(X1, y1, PM$M.lasso,    lambda,   alpha,   niter=20)

# Elastic Net Regularized Logistic Regression Model
out3 =    SGNLR(X1, y1, PM$M.elastic,  lambda,   alpha,   niter=20)

# Classical Network-regularized Logistic Regression Model
out4 =    SGNLR(X1, y1, PM$M.network,  lambda,   alpha,   niter=20)

# Adaptive Network-regularized Logistic Regression Model
out5 =    SGNLR(X1, y1, PM$M.AdaptNet, lambda,   alpha,   niter=20)

# Absolute Network-regularized Logistic Regression Model
out6 = abs.SNLR(X1, y1, PM$M.network,  lambda,   alpha,   niter=20)

# Testing
# X2 <- X1; y2 <- y1
res1 = predict.SGNLR(X2,y2,out1$w)
res2 = predict.SGNLR(X2,y2,out2$w)
res3 = predict.SGNLR(X2,y2,out3$w)
res4 = predict.SGNLR(X2,y2,out4$w)
res5 = predict.SGNLR(X2,y2,out5$w)
res6 = predict.SGNLR(X2,y2,out6$w)



# 保存 系数 -------------------------------------------------------------------

coef_lasso <- as.matrix(out2$w[-101])
rownames(coef_lasso) <- c(1:100)
sum(coef_lasso != 0)
# # View(coef_lasso)
# write.csv(coef_lasso, file = "theta_lasso.csv",row.names = F)

coef_elastic <- as.matrix(out3$w[-101])
rownames(coef_elastic) <- c(1:100)
sum(coef_elastic != 0)
# View(coef_elastic)
# write.csv(coef_elastic, file = "theta_elastic.csv",row.names = F)
 
coef_network <- as.matrix(out4$w[-101])
rownames(coef_network) <- c(1:100)
sum(coef_network != 0)
# View(coef_network)
# write.csv(coef_network, file = "theta_network.csv",row.names = F)

# 输入 CNet-RLR 系数 ----------------------------------------------------------

# CNet_RLR <- as.matrix(read.csv("theta1.csv", header = F))
# 
# # CNet_RLR[which(CNet_RLR[,1] <= 1e-8),] <- 0
# 
# library(pROC)
# X = X2;  y0 = y2;  w = CNet_RLR
# n = dim(X)[1]; p = dim(X)[2]
# if(!is.matrix(w)){w = matrix(c(w),ncol = 1)}
# # X = cbind(X,rep(1,n)) # X is a n x (p+1) matrix
# pp = 1/(1+exp(-X%*%w))
# y = rep(1,n)
# # View(cbind(y0,pp))
# y[pp<0.5] = 0
# accuracy.rate = length(which(y==y0))/n
# # AUC = accuracy.rate
# AUC = auc(y0, c(pp))[1]
# # AUC = roc(as.factor(y0), c(pp))$auc[1]
# accuracy.rate
# AUC

# 提取特征 --------------------------------------------------------------------

# coef_network <- as.matrix(out4$w[-101])
# rownames(coef_network) <- c(1:100)
# non_zero <- as.matrix(coef_network[-which(coef_network[,1] == 0),])
# feature <- rownames(non_zero)
# View(feature)
# write.csv(feature, file = "feature_network_example_lasso.csv",row.names = F)
# write.csv(feature, file = "feature_network_example_elastic.csv",row.names = F)
# write.csv(feature, file = "feature_network_example.csv",row.names = F)
## ===========================================================================

multi.method.AUC = c(res1$AUC,res2$AUC,res3$AUC,res4$AUC,res5$AUC,res6$AUC)
names(multi.method.AUC) = c("LR","Lasso.LR","Elastic.LR","Network.LR","AdaNet.LR","AbsNet.LR")

# Result in term of AUC
print(multi.method.AUC)

sum(coef_lasso != 0)
sum(coef_elastic != 0)
sum(coef_network != 0)
