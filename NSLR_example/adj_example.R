
# 试图将 模拟数据 中的邻接矩阵 转换成 连接对 -------------------------------------------------

rm(list=ls())

library(tidyfst)

# 读入数据 --------------------------------------------------------------------

setwd('C:/Users/LiLingyu/Desktop/NSLR-master')   
adj <- as.matrix(read.csv('adj_example.csv'))
label <- c(1:100)

colnames(adj) <- label
rownames(adj) <- label

adj[upper.tri(adj)] <- 0
adj_up <- adj

# 模拟数据 --------------------------------------------------------------------

tdf_list = mat_df(adj_up)
dim(tdf_list)    # adj_up
head(tdf_list)

list1 <- tdf_list[-which(tdf_list[,3] == "0"),]
list2 <- list1[,-3]
# write.csv(list2, file = "list.csv", row.names = F, quote = F)
