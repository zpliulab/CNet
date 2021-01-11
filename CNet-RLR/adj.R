library(igraph)

## 输入gene 和gene数对，得到邻接矩阵

rm(list=ls())

setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')
# gene <- read.csv('gene_net_neworder_cuttwo.csv')
gene <- read.csv('gene_net_neworder.csv')
genes <- data.frame(gene[,1])
# net <- read.csv('net_in_genes_adjp_cuttwo.csv')
net <- read.csv('net_in_genes_adjp.csv')

G <- graph_from_data_frame(net, directed=F, vertices=genes)
print(G, e=TRUE, v=TRUE)
# plot(G)


# 将图转换为邻接矩阵 ---------------------------------------------------------------

# adj <- as_adjcaency_matrix(G,sparse=FALSE)  # 作用同 get.adjacency
adj <- get.adjacency(G,sparse=FALSE) 
# write.csv(adj, 'adj_net.csv')
# write.csv(adj, 'adj_net_cuttwo.csv')



# 拉普拉斯矩阵 ------------------------------------------------------------------

# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj             # 最普通的 L 矩阵 
  return(L)
}

L <- Non.NormalizedLaplacianMatrix(adj)
# 特征值 ---------------------------------------------------------------------

a.e <- eigen(L,symmetric=T)

Vector <- a.e$vectors
eigvalue <- a.e$values

View(a.e$vectors)
a.e$vectors
View(a.e$values)
# setwd('D:\\E\\博士\\R_程序\\UCEC\\Data\\no_cut')
# write.csv(Vector, 'Vector_R.csv')
# write.csv(eigvalue, 'eigvalue_R.csv')
# write.csv(Vector, 'Vector_cuttwo_R.csv')
# write.csv(eigvalue, 'eigvalue_cuttwo_R.csv')

# 归一化的拉普拉斯矩阵 --------------------------------------------------------------

# Normalized Laplacian Matrix from adjacency matrix
laplacianMatrix = function(adj){
  diag(adj) <- 0                   # 邻接矩阵对角元0
  # 度矩阵元素（对角）--邻接矩阵每行元素的绝对值之和 
  deg <- apply(abs(adj),1,sum)     # abs(adj)-矩阵各元素去绝对值、1-表示按行计算，2表示按列、sum-自定义的调用函数
  p <- ncol(adj)
  L <- matrix(0,p,p)               # p*p 的0元素的 Laplaceian 矩阵
  nonzero <- which(deg!=0)         # 哪些行 元素绝对值之和不为0
  for (i in nonzero){
    for (j in nonzero){
      L[i,j] <- -adj[i,j]/sqrt(deg[i]*deg[j])  # i j 不等时（L 为对称阵）
    }
  }
  diag(L) <- 1                                 # 对角线为1
  return(L)
}

L_norm <- laplacianMatrix(adj)
# 特征值 ---------------------------------------------------------------------

a.e <- eigen(L_norm,symmetric=T)

Vector <- a.e$vectors
eigvalue <- a.e$values

View(a.e$vectors)
a.e$vectors
View(a.e$values)
# setwd('D:\\E\\博士\\R_程序\\UCEC\\Data\\no_cut')
# write.csv(Vector, 'Vector_R_norm.csv')
# write.csv(eigvalue, 'eigvalue_R_norm.csv')
# write.csv(Vector, 'Vector_cuttwo_R_norm.csv')
# write.csv(eigvalue, 'eigvalue_cuttwo_R_norm.csv')
