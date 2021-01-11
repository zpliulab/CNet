
###############   theta_selection (��matlab�����ϵ��theta��ѡ���������)  ##############

rm(list=ls())

setwd('D:\\E\\��ʿ\\R_����\\UCEC')               #�趨����Ŀ¼ΪD��
myfile = list.files("Data_theta126e50102")                #list.files���input�ļ����������ļ�������a
dir = paste("./Data_theta126e50102/", myfile, sep="")     #��paste�����·������dir
n = length(dir)                                  #��ȡdir���ȣ�Ҳ�����ļ����µ��ļ�����
theta = read.csv(file = dir[1],header=F,sep=",") #�����һ���ļ����ݣ����Բ����ȶ�һ��������Ϊ�˼򵥣�ʡȥ����data.frame��ʱ�䣬��ѡ���ȶ���һ���ļ���
for (i in 2:n){
  new_theta = read.csv(file = dir[i], header=F, sep=",")
  theta = cbind(theta,new_theta)
}                                                #ѭ���ӵڶ����ļ���ʼ���������ļ�������ϵ�merge.data������
setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data')  
gene <- read.csv('gene_net_neworder.csv', header = T, sep=',')
rownames(theta) <- gene[,1]


# ɸѡϵ��>ĳ������ ---------------------------------------------------------------

k = length(theta) 
label <- which(abs(theta[,12]) >= 1e-5)
# for (i in 2:k){
#   # i = 2
#   new_label = which(abs(theta[,i]) >= 1e-5)
#   print(paste("�ڣ�",i))
#   print(length(new_label))
#   # label = union(label,new_label)
#   label = intersect(label,new_label)
# }
feature <- gene[label,1]


# feature �������� ------------------------------------------------------------

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
plot(g, layout=layout.fruchterman.reingold, # ֻ����һ�У�ͼ������һ����
     vertex.size = 8,  # ���ýڵ��С
     vertex.label = V(g)$name, # ��Ȼ�ߺͽڵ���ܶ������֣���Ĭ��ʱ��Щ���ֿ���û�б�������ǩ
     vertex.label.cex = 0.8, # ��ǩ�����С
     vertex.label.dist = 0.1, # ���ýڵ�ͱ�ǩ�ľ��룬���ڴ����ص�
     vertex.label.color = "black"  # ���ñ�ǩ��ɫ
)
# write.csv(ed,"newresult0103\\net_in_feature_CNet-RLR.csv",row.names = F,quote = F)

# feature ɸѡ --------------------------------------------------------------

kk1 <- which( node_used %in% used[,1])    
kk2 <- which( node_used %in% used[,2])    
length(union(kk1,kk2))
feature_new <- as.matrix(node_used[union(kk1,kk2),])  # ϵ����Ϊ0�����������е�
# write.csv(feature_new,file = "newresult0103\\feature_CNet-RLR_0102.csv", row.names = F)


################ feature_selection (�����������fetature��Ӧ������) ###############
###################################################################################

library(dplyr)        
library(tidyr)
library(tidyverse)    

# gene -- feature; data -- ԭʼ�Ͷ�������
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

# ���� feature �� gene �������� --------------------------------------------------
Data1 = read.table("D:\\E\\��ʿ\\R_����\\UCEC\\GSE63678\\GSE63678_scale_ucec.txt", header = T, check.names = FALSE)
genedata2 <- my_feature_data(feature_new, Data1)
# write.table(genedata2,"newresult0103\\Feature63678_CNet-RLR_22_0102.txt",quote=F,sep="\t")


# ��ԭ���ݼ���ȡ���� ---------------------------------------------------------------
Data2 = read.table("UCEC_outcome_scale_net.txt", header = T, check.names = FALSE)
gene_overlap <- as.matrix(rownames(genedata2)[-1])
genedata_ucec <- my_feature_data(gene_overlap, Data2)
# write.table(genedata_ucec,"newresult0103\\Feature_ucec_Net-RLR_15_0102.txt",quote=F,sep="\t")
