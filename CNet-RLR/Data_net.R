

## 2020.7.16  ÿһ��������������ļ��Ľ��������ÿһ���ֶ������µ��õı����ļ���

setwd('D:\\E\\��ʿ\\R_����\\UCEC\\Data')


# ԭʼ����+��ǩ -----------------------------------------------------------------
data <- read.table('HiSeqV2.txt',header = T, check.names = FALSE)
lab <- read.table('label.txt',check.names = FALSE)
data_lab <- cbind(lab, t(data)) 
data_lab_t <- t(data_lab)
rownames(data_lab_t) <-data_lab_t[,1] 
Data <- data_lab_t[,-1]
# View(Data[1:10,1:10])
# write.table(Data, file = "UCEC_outcome.txt", quote=F, sep="\t")   


# �������� --------------------------------------------------------------------

Data <- read.table('UCEC_outcome.txt',header = T, check.names = FALSE)
# View(Data(1:10,1:10])
Data_label <- Data[1,]

# my_scale -------------------------------------------------------------------

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

# ���� ----------------------------------------------------------------------

Data1 <- my_scale(Data)
# mean(Data1[3,])     # ��ֵ�Ƿ�Ϊ0
# write.table(Data1,file = "UCEC_scale.txt",quote=F,sep="\t")  


# gene ��������  ��  net_sig ����  merge -----------------------------------------

Data1 <- read.table('UCEC_scale.txt',header = T, check.names = FALSE)
net_sig <- read.csv("renewed_Regnetwork_10_5_sig.csv",header = T)  # 20606


# �� Data1 ����������
Data2 <- cbind(rownames(Data1[-1,]), Data1[-1,])
colnames(Data2) <- c("Gene_symbol", colnames(Data1))
# View(Data2[1:10,1:10])

# �ϲ����� 
Data_net <- merge(net_sig, Data2, by.x="net_sig", by.y = "Gene_symbol", all=FALSE)
# dim(Data_net)    # 16245   202
rownames(Data_net) <- Data_net[,1]
Data_net_1 <- Data_net[,-1]
# View(Data_net_1[1:10,1:10])
# dim(Data_net_1)    # 16245   201

# ���ӱ�ǩ
Data_label <- Data1[1,]
Data_net_outcome <- rbind(Data_label, Data_net_1)
# View(Data_net_outcome[1:10,1:10])
# dim(Data_net_outcome)    # 16246   201
# write.table(Data_net_outcome,"UCEC_outcome_scale_net.txt",quote=F,sep="\t")  