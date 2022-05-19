
#5-24-2015
#prostate cancer

library(glasso)
t_static<- read.csv("t_static.csv",header=T)
t_static[,1]<-NULL
covA <- cov(t_static)
y <- seq(from=0.39, to=0.40, by=0.01) 
covR <- glassopath(covA,rholist=y) 
for (i in 1:2) {
  temp <- covR$wi[,,i] #wi:Estimated inverse covariance matrix, an array of dimension (nrow(s),ncol(n),length(rholist))
  write.table(x=temp, file=paste("t_static_matrix", i, sep = ""), sep="\t", row.names = F, col.names = F, quote = F)
}

genes<-colnames(t_static)

colnames(a1)<-genes
rownames(a1)<-genes
colnames(a2)<-genes
rownames(a2)<-genes

genepair(a1,"t_static_0.39.txt")
uniqgene("t_static_0.39.txt","uniq_t_static_0.39.txt")

genepair(a2,"t_static_0.40.txt")
uniqgene("t_static_0.40.txt","uniq_t_static_0.40.txt")

a39<-read.table("t_static_0.39.txt",header=F)
a39_col2<-as.data.frame(a39[,2])
colnames(a39_col2)<-c("Gene ID")
b39<-unique(a39_col2)
b39t<-t(b39)
write.table(b39t,"b39t.txt",col.names = F,row.names = T,sep=",",quote = T)

a40<-read.table("t_static_0.40.txt",header=F)
a40_col2<-as.data.frame(a40[,2])
colnames(a40_col2)<-c("Gene ID")
b40<-unique(a40_col2)
b40t<-t(b40)
write.table(b40t,"b40t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE7701_10 GSE4196911_650.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,1:10]
t_time<-data.frame(t(time))
time0.39<-t_time[,c("ALOX15","ALOX15B","ANPEP","CNTNAP2","ERG","F5","GDF15",
                    "GHR","HPGD","HPN","IL6","IL8","KRT5","PLA2G7","PLAUR",
                    "SERPINB5","TRIM29")]

time0.40<-t_time[,c("ALOX15","ALOX15B","ANPEP","ERG","F5","GHR","HPGD",
                    "IL6","IL8","KRT5","PLA2G7","PLAUR","SERPINB5",
                    "TRIM29")]
write.table(time0.40,"Time0.40_14genes.txt",col.names = T,row.names = F,sep="\t",quote = F)

#yeast
Total<-read.csv('GSE19213_4 GSE6212_11.csv',header=T)
rows<-Total[,1]
Total[,1]<-NULL
rownames(Total)<-rows

static<-Total[1:4]
t_static<-data.frame(t(static))
write.table(x =t_static, file = 'C:/Users/ZTQ/Desktop/t_static.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')

library(glasso)
t_static<- read.csv("t_static.csv",header=T)
t_static[,1]<-NULL
covA <- cov(t_static)
y <- seq(from=0.40, to=0.43, by=0.01) 
covR <- glassopath(covA,rholist=y) 
for (i in 1:4) {
  temp <- covR$wi[,,i] #wi:Estimated inverse covariance matrix, an array of dimension (nrow(s),ncol(n),length(rholist))
  write.table(x=temp, file=paste("t_static_matrix", i, sep = ""), sep="\t", row.names = F, col.names = F, quote = F)
}

# count non-zero number
library("Matrix")
a1<-read.table("t_static_matrix1",header=F) 
b1<-nnzero(a1,na.counted = TRUE)
print(b1-1266) 

a2<-read.table("t_static_matrix2",header=F) 
b2<-nnzero(a2,na.counted = TRUE)
print(b2-1266) 

a3<-read.table("t_static_matrix3",header=F) 
b3<-nnzero(a3,na.counted = TRUE)
print(b3-1266)

a4<-read.table("t_static_matrix4",header=F) 
b4<-nnzero(a4,na.counted = TRUE)
print(b4-1266)

genes<-rownames(static)
colnames(a3)<-genes
rownames(a3)<-genes

colnames(a2)<-genes
rownames(a2)<-genes

# unique gene pairs
genepair<-function(a,data_out){
  for (i in 1:nrow(a)) {
    col_index_t <- which(a[i,]!=0)
    length <- length(col_index_t)
    for (j in seq(from=1, to=length, by=1)){
      col_index<-col_index_t[j]
      col_name <- colnames(a)[col_index_t[j]]
      row_index <-i
      row_name <- rownames(a)[i]
      value <-a[i,col_index_t[j]]
      if (col_index!=row_index){
        data <- data.frame(row_index, row_name, col_index, col_name, value)
        write.table(data,data_out, sep="\t", row.names = F, col.names = F, quote = F, append = T)}
    }
  }
}
#Add column names and extract unique genes
uniqgene<-function(data_i,data_o)
{
  a<-read.table(data_i,header=F)
  name<-c("ID1","gene1","ID2","gene2","value")
  colnames(a)<-name
  a[,"value"]<-round(a[,"value"],digits =4 )
  b<-a[!duplicated(a[,c("value")]),]
  write.table(b,file =data_o,sep="\t", row.names = F, col.names = T, quote = F)
}


genepair(a3,"t_static_0.42.txt")
uniqgene("t_static_0.42.txt","uniq_t_static_0.42.txt")

genepair(a2,"t_static_0.41.txt")
uniqgene("t_static_0.41.txt","uniq_t_static_0.41.txt")

a42<-read.table("t_static_0.42.txt",header=F)
a42_col2<-as.data.frame(a42[,2])
colnames(a42_col2)<-c("Gene ID")
b42<-unique(a42_col2)
b42t<-t(b42)
write.table(b42t,"b42t.txt",col.names = F,row.names = T,sep=",",quote = T)

a41<-read.table("t_static_0.41.txt",header=F)
a41_col2<-as.data.frame(a41[,2])
colnames(a41_col2)<-c("Gene ID")
b41<-unique(a41_col2)
b41t<-t(b41)
write.table(b41t,"b41t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE19213_4 GSE6212_11.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,5:15]
t_time<-data.frame(t(time))
time0.41<-t_time[,c("ADA2","CGR1","DAK2","FMT1","GIT1","HOS2","MRPL4","RPA12",
                    "RRP12","SNU66","SPB1","SWF1","TEA1","VRP1","YOX1","ZDS1"
                    )]
write.table(time0.41,"Time0.41_16genes.txt",col.names = T,row.names = F,sep="\t",quote = F)
# 16 genes

time0.42<-t_time[,c("ADA2","CGR1","DAK2","FMT1","HOS2","RPA12","RRP12","SNU66",
                    "SPB1","SWF1","TEA1","VRP1","YOX1","ZDS1"
)]
write.table(time0.42,"Time0.42_14genes.txt",col.names = T,row.names = F,sep="\t",quote = F)
# 14 genes