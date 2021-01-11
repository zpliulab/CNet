
library(caret)
library(dplyr)

# ·�� ----------------------------------------------------------------------
setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data')

# ���� ----------------------------------------------------------------------
x <- read.table("UCEC_scale_net_adjp_net.txt", header = T, check.names = FALSE)
data <- data.frame(t(x))


# 0/1�ֿ� -------------------------------------------------------------------

data1 <- data[which(data$Label == 1),] 
data0 <- data[-which(data$Label == 1),] 

# ���� ���� ѵ����+���Լ� -----------------------------------------------------------

setwd('D:\\E\\��ʿ\\R_����\\UCEC')
dir.create("Data_train")
dir.create("Data_test")
for(i in 1: 20){
  # i = 1
  set.seed(i)

  n2<-nrow(data0) 
  data1_train <- data0[sample(n2,size = 24,replace = F), ] #���Ż������ȡ 24 ��������Ϊ��������
  # data1_test <-anti_join(data0,data1_train)                #A2���г�ȥtestdata�����ݵ�����
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
