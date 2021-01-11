############################# 下载 R 包 #############################
library(tidyverse)
library(caret)   # 十折交叉验证     
library(glmnet)
library(ggplot2)
library(caTools)
library(magrittr)
library(ROCR)
library(pROC)
library(glmnet)
library(ncpen)
library(stargazer) #转换latex 
library(broom) #回归结果保存
library(ncvreg)
library(l0ara)
library(stargazer)
library(plyr)
library(pROC)

rm(list=ls())

######################### glmnet 模拟 Lasso 惩罚 #############################

## 存储预测结果
coef_lasso <- matrix()       # 存储系数结果
pred_lasso <- matrix()   # 存储预测结果

## 训练集
setwd('D:\\E\\博士\\R_程序\\UCEC') 
myfile = list.files("Data_train")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_train/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
train <- as.matrix(read.table(file = dir[12],header=T))
x <- train[,-1]
y <- train[,1]

## 种子
i = 12
set.seed(i)

cv.lasso = cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10, type.measure = "class" )

## 拟合
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)
coef <- as.matrix(coef(lasso.model))
tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
step_lasso <- data.frame(cbind(coef, tcross))
coef_lasso <- cbind(coef_lasso, step_lasso)   #temp按行和pred合并

## 测试集
myfile = list.files("Data_test")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_test/", myfile, sep="")     #用paste命令构建路径变量dir
m = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
test <- as.matrix(read.table(file = dir[12],header=T))
x.test <- test[,-1]
y.test <- test[,1]

## 预测
p <- predict(lasso.model, newx = x.test, type = "response")
kcross <- rep(i, length(p)) 
temp_lasso <- data.frame(cbind(y.test, p, kcross)) # 真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
pred_lasso <- cbind(pred_lasso,temp_lasso)   #temp按行和pred合并

print(paste("第：",i)) 

# setwd("D:\\E\\博士\\R_程序\\UCEC\\Log_reg11\\lasso")
# write.csv(coef_lasso, file = "coef_lasso.csv")
# write.csv(pred_lasso, file = "pred_lasso.csv")


######################### glmnet 模拟 Elastic Net 惩罚 #############################

## 存储预测结果
coef_elastic <- matrix()   # 存储系数结果
pred_elastic <- matrix()   # 存储预测结果

## 训练集
setwd('D:\\E\\博士\\R_程序\\UCEC') 
myfile = list.files("Data_train")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_train/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
train <- as.matrix(read.table(file = dir[12],header=T))
x <- train[,-1]
y <- train[,1]

## 种子
i = 12
set.seed(i)

cv.elastic = cv.glmnet(x, y, alpha = 0.5, family = "binomial", nfolds = 10, type.measure = "class")

## 拟合
elastic.model <- glmnet(x, y, alpha = 0.5, family = "binomial", lambda = cv.elastic$lambda.1se)
coef <- as.matrix(coef(elastic.model))
tcross <- rep(i, length(coef))              # i是第几次循环交叉，共K次
step_elastic <- data.frame(cbind(coef, tcross))
coef_elastic <- cbind(coef_elastic, step_elastic)   #temp按行和pred合并

## 测试集
myfile = list.files("Data_test")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_test/", myfile, sep="")     #用paste命令构建路径变量dir
m = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
test <- as.matrix(read.table(file = dir[12],header=T))
x.test <- test[,-1]
y.test <- test[,1]

## 预测
p <- predict(elastic.model, newx = x.test, type = "response")
kcross <- rep(i, length(p)) 
temp_elastic <- data.frame(cbind(y.test, p, kcross))

pred_elastic <- cbind(pred_elastic,temp_elastic)   #temp按行和pred合并

print(paste("第：",i)) 

setwd('D:\\E\\博士\\R_程序\\UCEC\\Log_reg11\\elastic_net')
# write.csv(coef_elastic, file = "coef_elastic.csv")
# write.csv(pred_elastic, file = "pred_elastic.csv")

########################################  基因合并 ############################ 
############################################################################### 

rm(list=ls())
setwd("D:\\E\\博士\\R_程序\\UCEC\\Log_reg11")
Coef_lasso <- read.table("lasso\\coef_lasso.csv", header=TRUE, sep = ',')[,c(1,3)]
Coef_Elastic <- read.table("elastic_net\\coef_elastic.csv", header=TRUE, sep = ',')[,c(1,3)]

# 计算个数 --------------------------------------------------------------------
sum(Coef_lasso[-1,2] != 0)
sum(Coef_Elastic[-1,2] != 0)

lasso <- Coef_lasso[which(Coef_lasso[-1,2] != 0),1]
Elastic <- Coef_Elastic[which(Coef_Elastic[-1,2] != 0),1]

setwd('D:\\E\\博士\\R_程序\\UCEC\\Log_reg11\\result')
# write.csv(lasso, "coef_lasso1.csv", row.names = F)
# write.csv(Elastic, "coef_Elastic1.csv", row.names = F)

