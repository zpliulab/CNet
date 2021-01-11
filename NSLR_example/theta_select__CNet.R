
rm(list=ls())

setwd('C:\\Users\\LiLingyu\\Desktop\\NSLR-master')               #设定工作目录为D盘
myfile = list.files("Data_theta")                #list.files命令将input文件夹下所有文件名输入a
dir = paste("./Data_theta/", myfile, sep="")     #用paste命令构建路径变量dir
n = length(dir)                                  #读取dir长度，也就是文件夹下的文件个数
theta = read.csv(file = "theta1.csv",header=F,sep=",") #读入第一个文件内容（可以不用先读一个，但是为了简单，省去定义data.frame的时间，我选择先读入一个文件。


gene <- as.matrix(c(1:100))
rownames(theta) <- gene[,1]

# 筛选系数>某个数的 ---------------------------------------------------------------

k = length(theta) 
label <- which(abs(theta[,1]) >= 1e-5)
feature <- gene[label,1]

# write.csv(feature,file = "feature_CNet-RLR.csv", row.names = F)


# feature 在网络中 ------------------------------------------------------------

net_used <- as.matrix(read.csv("list.csv",header = T))
node_used <- as.matrix(feature)

k1 <- which(net_used[,1] %in% node_used)   # 160
k2 <- which(net_used[,2] %in% node_used)   #60
length(intersect(k1,k2))    # 115  13
used <- net_used[intersect(k1,k2),]

library(igraph)
PP <- graph_from_data_frame(used,directed = F)
p1 <- simplify(PP,remove.multiple = TRUE, remove.loops = TRUE)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)

g <- p1
plot(g, layout=layout.fruchterman.reingold, # 只有这一行，图都挤到一块了
     vertex.size = 8,  # 设置节点大小
     vertex.label = V(g)$name, # 虽然边和节点可能都有名字，但默认时这些名字可能没有被当做标签
     vertex.label.cex = 0.8, # 标签字体大小
     vertex.label.dist = 0.1, # 设置节点和标签的距离，便于错开重叠
     vertex.label.color = "black"  # 设置标签颜色
)
# write.csv(ed,"net_in_feature_CNEt-RLR.csv",row.names = F,quote = F)


# feature 筛选 --------------------------------------------------------------

kk1 <- which( node_used %in% used[,1])   # 5
kk2 <- which( node_used %in% used[,2])   # 28
length(union(kk1,kk2))
feature_new <- as.matrix(node_used[union(kk1,kk2),])
# write.csv(feature_new,file = "feature_in_net_CNet-RLR.csv", row.names = F)

