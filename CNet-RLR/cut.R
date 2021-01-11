# install.packages("igraph")
library(igraph)

setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')

gene <- read.csv('gene_net_neworder.csv')
genes <- data.frame(gene[,1])
# net <- read.csv('net_in_genes_adjp.csv')
# g <- graph_from_data_frame(net, directed=F, vertices=genes)

net <- read.csv('net_in_genes_adjp.csv')
g <- graph_from_data_frame(net, directed=F)

print(g, e=TRUE, v=TRUE)
plot(g, layout=layout.fruchterman.reingold, # 只有这一行，图都挤到一块了
     vertex.size=4,  # 设置节点大小
     vertex.label = V(g)$name, # 虽然边和节点可能都有名字，但默认时这些名字可能没有被当做标签
     vertex.label.cex=0.7, # 标签字体大小
     vertex.label.dist=0.4, # 设置节点和标签的距离，便于错开重叠
     vertex.label.color = "black"  # 设置标签颜色
)

# 计算距离最大的两个点 ------------------------------------------------------------------

## 直径，breadth-first search
diameter(g) 

## TRUE, the diameters of the connected components
diameter(g, unconnected=TRUE)

## FALSE, the number of vertices
diameter(g, unconnected=FALSE)

## Weighted diameter
# E(g)$weight <- sample(seq_len(ecount(g)))  # ecount 计算g的边数

## returns a path with the actual diameter
get_diameter(g) 
# path <- get_diameter(g) 
# plot(path)
# get_diameter(g, weights=NA)

## returns two vertex ids, connected by the diameter path.
farthest_vertices(g) 
# diameter(g, weights=NA)

# 找最小 cut  ----------------------------------------------------------------

## FALSE, the edges in the cut and a the two (or more) partitions are also returned.
min_cut(g, source = "ESPL1", target = "ZBTB16", value.only = FALSE)
# min_cut(g, source = 2, target = 5, value.only = FALSE)

## TRUE, only the minumum cut value is returned
min_cut(g, source = "ESPL1", target = "ZBTB16", value.only = TRUE)
# min_cut(g, source = 2, target = 5, value.only = TRUE)

# 组装矩阵和向量 -----------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')

library(igraph)
gene <- read.csv('gene_net_neworder.csv',header = T)

net <- read.csv('net_in_genes_adjp.csv')
g <- graph_from_data_frame(net, directed=F)

node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")


# 2020.7.20 输出1个向量，向量中定点为 1，割点位为 -1 ----------------------------------------------------

my_vector <- function(gene, lab1){
  k <- length(as.matrix(gene))
  vec1 <- matrix(0,1,k)
  
  l <- length(lab1)

  for (i in 1:l) {
    # i= 1
    vec1[which(gene[,1]==lab1[i])] <- c("-1")
  }
  
  x_first <- which(gene[,1]==lab1[1])
  x_end <- which(gene[,1]==lab1[l])
  vec1[,c(x_first,x_end)] <- c("1")
  
  return(vec1)
}

vec2 <- my_vector(gene, lab1)
vec3 <- t(vec2)
View(vec3)
# setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')
# write.table(vec3,file = "vector_hat1_1225.txt", quote=F, sep="\t", row.names = F)

