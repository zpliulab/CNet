library(igraph)

## ����gene ��gene���ԣ��õ��ڽӾ���

rm(list=ls())

setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data')
# gene <- read.csv('gene_net_neworder_cuttwo.csv')
gene <- read.csv('gene_net_neworder.csv')
genes <- data.frame(gene[,1])
# net <- read.csv('net_in_genes_adjp_cuttwo.csv')
net <- read.csv('net_in_genes_adjp.csv')

G <- graph_from_data_frame(net, directed=F, vertices=genes)
print(G, e=TRUE, v=TRUE)
# plot(G)


# ��ͼת��Ϊ�ڽӾ��� ---------------------------------------------------------------

# adj <- as_adjcaency_matrix(G,sparse=FALSE)  # ����ͬ get.adjacency
adj <- get.adjacency(G,sparse=FALSE) 
# write.csv(adj, 'adj_net.csv')
# write.csv(adj, 'adj_net_cuttwo.csv')



# ������˹���� ------------------------------------------------------------------

# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj             # ����ͨ�� L ���� 
  return(L)
}

L <- Non.NormalizedLaplacianMatrix(adj)
# ����ֵ ---------------------------------------------------------------------

a.e <- eigen(L,symmetric=T)

Vector <- a.e$vectors
eigvalue <- a.e$values

View(a.e$vectors)
a.e$vectors
View(a.e$values)
# setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data\\no_cut')
# write.csv(Vector, 'Vector_R.csv')
# write.csv(eigvalue, 'eigvalue_R.csv')
# write.csv(Vector, 'Vector_cuttwo_R.csv')
# write.csv(eigvalue, 'eigvalue_cuttwo_R.csv')

# ��һ����������˹���� --------------------------------------------------------------

# Normalized Laplacian Matrix from adjacency matrix
laplacianMatrix = function(adj){
  diag(adj) <- 0                   # �ڽӾ���Խ�Ԫ0
  # �Ⱦ���Ԫ�أ��Խǣ�--�ڽӾ���ÿ��Ԫ�صľ���ֵ֮�� 
  deg <- apply(abs(adj),1,sum)     # abs(adj)-�����Ԫ��ȥ����ֵ��1-��ʾ���м��㣬2��ʾ���С�sum-�Զ���ĵ��ú���
  p <- ncol(adj)
  L <- matrix(0,p,p)               # p*p ��0Ԫ�ص� Laplaceian ����
  nonzero <- which(deg!=0)         # ��Щ�� Ԫ�ؾ���ֵ֮�Ͳ�Ϊ0
  for (i in nonzero){
    for (j in nonzero){
      L[i,j] <- -adj[i,j]/sqrt(deg[i]*deg[j])  # i j ����ʱ��L Ϊ�Գ���
    }
  }
  diag(L) <- 1                                 # �Խ���Ϊ1
  return(L)
}

L_norm <- laplacianMatrix(adj)
# ����ֵ ---------------------------------------------------------------------

a.e <- eigen(L_norm,symmetric=T)

Vector <- a.e$vectors
eigvalue <- a.e$values

View(a.e$vectors)
a.e$vectors
View(a.e$values)
# setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data\\no_cut')
# write.csv(Vector, 'Vector_R_norm.csv')
# write.csv(eigvalue, 'eigvalue_R_norm.csv')
# write.csv(Vector, 'Vector_cuttwo_R_norm.csv')
# write.csv(eigvalue, 'eigvalue_cuttwo_R_norm.csv')