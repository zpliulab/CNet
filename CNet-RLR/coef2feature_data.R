
###############   theta_selection (从matlab输出的系数theta中选择出特征来)  ##############

rm(list=ls())

setwd('D:\\E\\博士\\R_程序\\UCEC')               #设定工作目录为D盘
myfile = list.files("Data_theta126e50102")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_theta126e50102/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
theta = read.csv(file = dir[1],header=F,sep=",") #读入第一个文件内容（可以不用先读一个，但是为了简单，省去定义data.frame的时间，我选择先读入一个文件。
for (i in 2:n){
  new_theta = read.csv(file = dir[i], header=F, sep=",")
  theta = cbind(theta,new_theta)
}                                                #循环从第二个文件开始读入所有文件，并组合到merge.data变量中
setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')  
gene <- read.csv('gene_net_neworder.csv', header = T, sep=',')
rownames(theta) <- gene[,1]


# 筛选系数>某个数的 ---------------------------------------------------------------

k = length(theta) 
label <- which(abs(theta[,12]) >= 1e-5)
# for (i in 2:k){
#   # i = 2
#   new_label = which(abs(theta[,i]) >= 1e-5)
#   print(paste("第：",i))
#   print(length(new_label))
#   # label = union(label,new_label)
#   label = intersect(label,new_label)
# }
feature <- gene[label,1]


# feature 在网络中 ------------------------------------------------------------

net_used <- as.matrix(read.csv("net_in_genes_adjp.csv",header = T))
node_used <- as.matrix(feature)

k1 <- which(net_used[,1] %in% node_used) 
k2 <- which(net_used[,2] %in% node_used) 
length(intersect(k1,k2))   
used <- net_used[intersect(k1,k2),]

library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP)
# p1 <- simplify(PP,remove.multiple = TRUE, remove.loops = TRUE) 
ed <- as_edgelist(p1, names = TRUE)

g <- p1
plot(g, layout=layout.fruchterman.reingold, # 只有这一行，图都挤到一块了
     vertex.size = 8,  # 设置节点大小
     vertex.label = V(g)$name, # 虽然边和节点可能都有名字，但默认时这些名字可能没有被当做标签
     vertex.label.cex = 0.8, # 标签字体大小
     vertex.label.dist = 0.1, # 设置节点和标签的距离，便于错开重叠
     vertex.label.color = "black"  # 设置标签颜色
)
# write.csv(ed,"newresult0103\\net_in_feature_CNet-RLR.csv",row.names = F,quote = F)

# feature 筛选 --------------------------------------------------------------

kk1 <- which( node_used %in% used[,1])    
kk2 <- which( node_used %in% used[,2])    
length(union(kk1,kk2))
feature_new <- as.matrix(node_used[union(kk1,kk2),])  # 系数不为0，且在网络中的
# write.csv(feature_new,file = "newresult0103\\feature_CNet-RLR_0102.csv", row.names = F)


################ feature_selection (从数据中输出fetature对应的数据) ###############
###################################################################################

library(dplyr)        
library(tidyr)
library(tidyverse)    

# gene -- feature; data -- 原始和独立数据
my_feature_data <- function(gene, Data1){   
  colnames(gene) <- c('gene')
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  
  genedata <- merge(gene, Data2, by = "gene")
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[1,],genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  
  return(genedata2)
}

# 根据 feature 找 gene 表达数据 --------------------------------------------------
Data1 = read.table("D:\\E\\博士\\R_程序\\UCEC\\GSE63678\\GSE63678_scale_ucec.txt", header = T, check.names = FALSE)
genedata2 <- my_feature_data(feature_new, Data1)
# write.table(genedata2,"newresult0103\\Feature63678_CNet-RLR_22_0102.txt",quote=F,sep="\t")


# 在原数据集提取特征 ---------------------------------------------------------------
Data2 = read.table("UCEC_outcome_scale_net.txt", header = T, check.names = FALSE)
gene_overlap <- as.matrix(rownames(genedata2)[-1])
genedata_ucec <- my_feature_data(gene_overlap, Data2)
# write.table(genedata_ucec,"newresult0103\\Feature_ucec_Net-RLR_15_0102.txt",quote=F,sep="\t")

