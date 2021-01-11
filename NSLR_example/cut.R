# install.packages("igraph")
library(igraph)

setwd('C:/Users/LiLingyu/Desktop/NSLR-master') 
gene <- as.matrix(c(1:100))
genes <- data.frame(gene[,1])

net <- read.csv('list.csv')
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

## returns a path with the actual diameter
get_diameter(g) 

## returns two vertex ids, connected by the diameter path.
farthest_vertices(g) 


# 找最小 cut  ----------------------------------------------------------------

## TRUE, only the minumum cut value is returned
min_cut(g, source = "36", target = "3", value.only = TRUE)
# min_cut(g, source = 2, target = 5, value.only = TRUE)

# 组装矩阵和向量 -----------------------------------------------------------------

library(igraph)
g <- graph_from_data_frame(net, directed=F)

node <- as.matrix(get_diameter(g))
lab1 <- as.matrix(rownames(node))
colnames(lab1) <- c("node")


# 2020.7.20 输出1个向量，向量为定点和割点位的系数 ----------------------------------------------------

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

# write.table(vec3,file = "vector_hat_example.txt", quote=F, sep="\t", row.names = F)

