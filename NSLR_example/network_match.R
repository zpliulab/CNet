## 将gene输进去，得到pathway gene数对

rm(list=ls())
setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')
net <- as.matrix(read.csv("renewed_Regnetwork_10_5.csv",header = T))
net[1,]

# node <- as.matrix(read.csv("gene_net_neworder_cuttwo.csv",header = T))#[-1,]
# node <- as.matrix(read.csv("gene_net_neworder.csv",header = T))#[-1,]

## 2020.1.2
node <- as.matrix(read.csv("gene_net_neworder_0102.csv",header = T))#[-1,]


node[1]
node_used <- node
net_used <- net[,c(2,4)]
k1 <- which(net_used[,1] %in% node_used)   # 4701    2042
k2 <- which(net_used[,2] %in% node_used)   # 7133     957
length(intersect(k1,k2))    # 230  13      # 193
# length(union(k1,k2))
used <- net_used[intersect(k1,k2),]    # 192
# used <- net_used[union(k1,k2),]

library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP,remove.loops = T,remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)
# write.csv(ed,"net_in_genes_adjp_cuttwo.csv",row.names = F,quote = F)


# 计算图的最大(弱或强)连通分量 ---------------------------------------------------------

# g <- sample_gnp(20, 1/20)
clu <- components(p1)
groups(clu)

g <- p1

# pdf(file = "net_in_genes_adjp_cuttwo.pdf",width = 15,height = 15)
plot(g, layout=layout.fruchterman.reingold, # 只有这一行，图都挤到一块了
     vertex.size=4,  # 设置节点大小
     vertex.label = V(g)$name, # 虽然边和节点可能都有名字，但默认时这些名字可能没有被当做标签
     vertex.label.cex=0.7, # 标签字体大小
     vertex.label.dist=0.4, # 设置节点和标签的距离，便于错开重叠
     vertex.label.color = "black"  # 设置标签颜色
     )
# dev.off()

eb <- cluster_edge_betweenness(g)
eb[[1]]
eb[[2]]
eb[[3]]
