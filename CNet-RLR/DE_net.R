## limma 分析得到差异gene

# install.packages("modelr")
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用
library(limma)
library(edgeR)

setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')

Data = read.table("UCEC_outcome_scale_net.txt", header = T, check.names = FALSE)
# View(Data[1:10,1:10])


# 表型标签 ----------------------------------------------------------------------

sample <- as.matrix(colnames(Data))
label <- t(Data[1,])
class <- cbind(sample, label)
colnames(class) <- c("sample","class")
# write.csv(class, file = "phenotype.csv", row.names = F)


# 排序 后 重组数据----------------------------------------------------------------------

Data1 <- Data[,which(Data[1,] == 1)]
# View(Data1[1:10,1:10])
Data0 <- Data[,which(Data[1,] == 0)]
# View(Data0[1:10,1:10])
data <- cbind(Data1,Data0)[-1,]
# View(data[1:10,1:10])

phe = read.csv("phenotype.csv", header=TRUE, sep = ',')
phe1 <- phe[which(phe[,2] == 1),]
phe1[,2] <- c("Normal")
phe0 <- phe[which(phe[,2] == 0),]
phe0[,2] <- c("Tumor")
phenotype <- rbind(phe1,phe0)
# write.table(phenotype, file = "phe.txt", quote=F, sep="\t", row.names = F)


# 差异分析 ----------------------------------------------------------------
disease = read.table("phe.txt",header = TRUE, sep = "\t")
disease <- factor(disease[,"class"])

# 构建实验设计矩阵
design <- model.matrix(~-1+disease)
# 构建对比模型,比较两个实验条件下表达数据
contrast.matrix <- makeContrasts(contrasts = "diseaseNormal - diseaseTumor", levels = design)
# 线性模型拟合
fit <- lmFit(data,design)
# 根据对比模型进行差值计算
fit1 <- contrasts.fit(fit,contrast.matrix)
# 贝叶斯检验
fit2 <- eBayes(fit1)
# 生成所有基因的检验结果报表
dif <- topTable(fit2,coef = "diseaseNormal - diseaseTumor",n = nrow(fit2),lfc = log2(2))
# dif <- topTable(fit2,coef = "diseaseNormal - diseaseTumor",n = nrow(fit2))
# write.csv(dif, file = "dif_FC_4324.csv")
# dif <- dif[dif[,"P.Value"]<0.05,]
dif0.05 <- dif[dif[,"adj.P.Val"] < 0.05,]    # 4324
dif_FC <- dif0.05[abs(dif0.05[,"logFC"]) > 2,]     # 360  # 1553
# View(dif0.05)
# write.csv(dif_FC, file = "dif_FC.csv")

gene_adjp <- rownames(dif_FC)
# write.csv(gene_adjp, file = "gene_adjp.csv", row.names = F)


# 挑选只在网络中的 gene  的表达数据 ----------------------------------------------------


## 2020.1.2   
net <- as.matrix(read.csv("renewed_Regnetwork_10_5.csv",header = T))
net_used <- net[,c(2,4)]

# net_used <- as.matrix(read.csv("net_in_genes_adjp.csv",header = T))
node_used <- as.matrix(read.csv("gene_adjp.csv",header = T))#[-1,]

k1 <- which(net_used[,1] %in% node_used)   # 178
k2 <- which(net_used[,2] %in% node_used)   # 178
length(intersect(k1,k2))    # 178
used <- net_used[intersect(k1,k2),]

kk1 <- which( node_used %in% used[,1])   # 14
kk2 <- which( node_used %in% used[,2])   # 121
length(union(kk1,kk2))      # 126
feature <- as.matrix(node_used[union(kk1,kk2),])
# write.csv(feature,file = "gene_net.csv", row.names = F)
# write.csv(feature,file = "gene_net_0102.csv", row.names = F) 



# 差异gene的表达数据 -------------------------------------------------------------

# merge
data <- read.table("UCEC_outcome_scale_net.txt", header = T, check.names = FALSE)
data1 <- cbind(rownames(data),data)
colnames(data1) <- c("gene_symbol", colnames(data))
data2 <- data1[-1,]    # 删除第一行label

gene_adjp = read.csv("gene_net_0102.csv", header=TRUE, sep = ',')
colnames(gene_adjp) <- c("gene_adjp")

data3 <- merge(gene_adjp, data2, by.x="gene_adjp",by.y = "gene_symbol",all=FALSE) 
data4 <- rbind(data1[1,-1],data3[,-1])
row.names(data4) <- c("Label",as.character(data3[,1]))
# View(data4[1:10,1:10])
dim(data4)    # 361 201  # 127 201
# write.table(data4, file = "UCEC_scale_net_adjp_net.txt",quote = F, sep = "\t")


# 保存新的gene列表，与数据集对应的，然后去找网络 -----------------------------------------------

gene_adjp_neworder <- data3[,1]
View(gene_adjp_neworder)
# write.csv(gene_adjp_neworder, file = "gene_net_neworder.csv", row.names = F)
# write.csv(gene_adjp_neworder, file = "gene_net_neworder_0102.csv", row.names = F)
