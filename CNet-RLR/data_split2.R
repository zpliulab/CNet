
library(caret)
library(dplyr)

# 路径 ----------------------------------------------------------------------
setwd('D:\\E\\博士\\R_程序\\UCEC\\Data')

# 数据 ----------------------------------------------------------------------
x <- read.table("UCEC_scale_net_adjp_net.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))


# 0/1分开 -------------------------------------------------------------------

data1 <- data[which(data$Label == 1),] 
data0 <- data[-which(data$Label == 1),] 

# 数据 —— 训练集+测试集 -----------------------------------------------------------

setwd('D:\\E\\博士\\R_程序\\UCEC')
dir.create("Data_train")
dir.create("Data_test")
for(i in 1: 20){
  # i = 1
  set.seed(i)

  n2<-nrow(data0) 
  data1_train <- data0[sample(n2,size = 24,replace = F), ] #不放回随机抽取 24 个样本作为测试样本
  # data1_test <-anti_join(data0,data1_train)                #A2表中除去testdata中数据的样本
  data_new <- rbind(data1_train, data1)
  
  training.samples <- data_new$Label %>% createDataPartition(p = 0.75, list = FALSE)
  train.data  <- data_new[training.samples, ]
  test.data <- data_new[-training.samples, ]
  
  name <- as.character(i)
  train <- as.matrix(train.data)
  test <- as.matrix(test.data)
  
  path <- paste("./Data_train/",paste(name,".txt",sep=""),sep="")
  write.table(train,path,quote = F)
  
  path <- paste("./Data_test/",paste(name,".txt",sep=""),sep="")
  write.table(test,path,quote = F)
}

