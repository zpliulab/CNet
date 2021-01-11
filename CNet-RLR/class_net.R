library(pROC)

# �������ݼ� -------------------------------------------------------------------

rm(list=ls())

setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data\\newresult0103')

# Net-RLR
Data  = read.table("Feature63678_Net-RLR_15_0102.txt", header = T, check.names = FALSE)
Data2 = read.table("Feature_ucec_Net-RLR_15_0102.txt", header = T, check.names = FALSE)

# CNet-RLR
# Data  = read.table("Feature63678_CNet-RLR_22_0102.txt", header = T, check.names = FALSE)
# Data2 = read.table("Feature_ucec_CNet-RLR_22_0102.txt", header = T, check.names = FALSE)

# ## lasso
# Data  = read.table("Feature63678_lasso_5_0102.txt", header = T, check.names = FALSE)
# Data2 = read.table("Feature_ucec_lasso_5_0102.txt", header = T, check.names = FALSE)

# ## elastic
# Data  = read.table("Feature63678_elastic_24_0102.txt", header = T, check.names = FALSE)
# Data2 = read.table("Feature_ucec_elastic_24_0102.txt", header = T, check.names = FALSE)

###################################################################################

x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]

# glm.fit <- glm(y.train~., data = x.train, family = binomial)
glm.fit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
summary(glm.fit)

# train --------------------------------------------------------------------

p_train <- predict(glm.fit, x.train, type = "response")
pred_glm <- cbind(p_train, y.train)
colnames(pred_glm) <- c('y.train', 'Lable')

p <- pred_glm[,1]
p_glm <- cbind(log(p/(1-p)), pred_glm[,2])
colnames(p_glm) <- c('y.train', 'Lable')

library(pROC)
p_train = as.matrix(p_train)
A_train <- data.frame(p_train, y.train)
names(A_train)<- c("p", "outcome")

# pdf(file = "pAUC_test_glm_cv4.pdf")
plot.roc(A_train$outcome, A_train$p, print.auc=T, main="pAUC")
# plot.roc(A_test$outcome, A_test$p)
# legend("bottomright", legend=c("Acc=0.875", "Pre=1.000 ", "Sn=0.75", "F-measure=0.857", "Sp=1.000", "AUC=0.875"))
# dev.off()

# test --------------------------------------------------------------------

p_test <- predict(glm.fit, x.test, type = "response")
pred_glm <- cbind(p_test, y.test)
colnames(pred_glm) <- c('y.test', 'Lable')

p <- pred_glm[,1]
p_glm <- cbind(log(p/(1-p)), pred_glm[,2])
colnames(p_glm) <- c('y.test', 'Lable')
# write.table(p_glm,"pre63678log.txt",quote=F,sep="\t")

library(pROC)
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
names(A_test)<- c("p", "outcome")

# pdf(file = "pAUC_test_63678.pdf",width = 5,height = 5)
plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
# plot.roc(A_test$outcome, A_test$p)
# legend("bottomright", legend=c("Acc = 0.833", "Pre = 0.800 ", "Sn = 0.800", 
#                                "Sp = 0.857", "F-measure = 0.800", "AUC = 0.914"))
# dev.off()


# ����4��RLR��AUC���� -----------------------------------------------------------

# ROC <- A_test
# ROC <- cbind(ROC, A_test)
# ROC1 <- ROC[,c(2,1,3,5,7)]
# write.csv(ROC1, file = "ROC_RLR.csv", row.names = F)


## ����

pred_glm <- A_test
predict = ifelse(pred_glm[,1] > 0.5, 1, 0)
predict_value = predict
true_value = pred_glm[,2]
error = predict_value-true_value

data <- t(Data)
# ����ģ��׼ȷ�ԣ�accuracy������ȷ�ȣ�Precision�����ٻ��ʣ�Recall��-- �����ԣ�sensitivity��-- �������ʣ�TPR���� F��ȣ�F-measure���ͻ�������
accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
precision = sum(true_value & predict_value)/sum(predict_value)  #��ʵֵԤ��ֵȫΪ1 / Ԥ��ֵȫΪ1 --- ��ȡ������ȷ��Ϣ����/��ȡ������Ϣ����
recall = sum(predict_value & true_value)/sum(true_value)        #��ʵֵԤ��ֵȫΪ1 / ��ʵֵȫΪ1 --- ��ȡ������ȷ��Ϣ���� /�����е���Ϣ����
# P��Rָ����ʱ�����ֵ�ì�ܵ��������������Ҫ�ۺϿ������ǣ�����ķ�������F-Measure���ֳ�ΪF-Score��
F_measure = 2*precision*recall/(precision+recall)    #F-Measure��Precision��Recall��Ȩ����ƽ������һ���ۺ�����ָ��
specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))

accuracy
precision
recall
specificity
F_measure
## ����������ʾ�������ΪTP��FN��FP��TN
table(true_value, predict_value) 