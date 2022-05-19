library(Biobase)
library(GEOquery)
gset<-getGEO("GSE770", GSEMatrix =TRUE)
gset<-gset[[1]]
show(gset)
str(gset)
e<-data.frame(exprs(gset))
boxplot(log(e))
# Already normolized

symb<-gset@featureData @data $"Gene Symbol"  
symb<-as.character(symb)
# We don't actually have gene symbols for all rows, we need to address for it
sel<-apply(as.matrix(symb,ncol=1),1,
           function(x){if(x=="")return(FALSE)else return(TRUE)})
levels(as.factor(sel))
selexpr<-c()
selsymb<-c()
for (i in 1:nrow(e)){
  if (sel[i]){
    selexpr<-rbind(selexpr,e[i,])
    selsymb<-rbind(selsymb,symb[i])
  }
}
selexpr<-data.frame(selexpr)
# get average value for each gene
ave<-aggregate(selexpr, by=list(selsymb), 
               FUN=mean)
rows<-ave[,1]
ave[,1]<-NULL
rownames(ave)<-rows
write.table(x =ave, file ="ave_GSE770.csv", sep =',', col.names=NA,row.names=T, quote = F, na = '')
# get the same genes from time-series data and static data
static<-read.csv('ave_allcancer_Static.csv',header=T)
m1 <- merge(ave, static, by.x = "row.names", by.y = "X")
rows<-m1[,1]
m1[,1]<-NULL
rownames(m1)<-rows
write.table(x =m1, file = 'GSE7701_10 GSE4196911_650.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')




m<-read.csv('GSE7701_10 GSE4196911_650.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,1:10]
static<-m[,11:649]


# log2 transform
qx <- as.numeric(quantile(static, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { static[which(static <= 0)] <- NaN
            static_log <- log2( static) }
t_static<-data.frame(t(static_log))
write.table(x =t_static, file = 'C:/Users/ZTQ/Desktop/t_static.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')

library(glasso)
t_static<- read.csv("t_static.csv",header=T)
t_static[,1]<-NULL
covA <- cov(t_static)
y <- seq(from=0.28, to=0.34, by=0.02) 
covR <- glassopath(covA,rholist=y) 
for (i in 1:4) {
  temp <- covR$wi[,,i] #wi:Estimated inverse covariance matrix, an array of dimension (nrow(s),ncol(n),length(rholist))
  write.table(x=temp, file=paste("t_static_matrix", i, sep = ""), sep="\t", row.names = F, col.names = F, quote = F)
}

# count non-zero number
library("Matrix")
a1<-read.table("t_static_matrix1",header=F) 
b1<-nnzero(a1,na.counted = TRUE)
print(b1-424) #0.28  128

a2<-read.table("t_static_matrix2",header=F) 
b2<-nnzero(a2,na.counted = TRUE)
print(b2-424) #0.3 92

a3<-read.table("t_static_matrix3",header=F) 
b3<-nnzero(a3,na.counted = TRUE)
print(b3-424) #0.32 72

a4<-read.table("t_static_matrix4",header=F) 
b4<-nnzero(a4,na.counted = TRUE)
print(b4-424) #0.34 46

genes<-rownames(static)
colnames(a1)<-genes
rownames(a1)<-genes

colnames(a2)<-genes
rownames(a2)<-genes

colnames(a3)<-genes
rownames(a3)<-genes

colnames(a4)<-genes
rownames(a4)<-genes

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
  a[,"value"]<-round(a[,"value"],digits =5 )
  b<-a[!duplicated(a[,c("value")]),]
  write.table(b,file =data_o,sep="\t", row.names = F, col.names = T, quote = F)
}

genepair(a1,"t_static_0.28.txt")
uniqgene("t_static_0.28.txt","uniq_t_static_0.28.txt")

genepair(a2,"t_static_0.3.txt")
uniqgene("t_static_0.3.txt","uniq_t_static_0.3.txt")

genepair(a3,"t_static_0.32.txt")
uniqgene("t_static_0.32.txt","uniq_t_static_0.32.txt")

genepair(a4,"t_static_0.34.txt")
uniqgene("t_static_0.34.txt","uniq_t_static_0.34.txt")

#abstract unique gene names from the files

a28<-read.table("t_static_0.28.txt",header=F)
a28_col2<-as.data.frame(a28[,2])
colnames(a28_col2)<-c("Gene ID")
b28<-unique(a28_col2)
b28t<-t(b28)
write.table(b28t,"b28t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE7701_10 GSE4196911_650.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,1:10]
t_time<-data.frame(t(time))
time0.28<-t_time[,c("ABCA5","ADIPOQ","ALOX12B","ALOX15","ALOX15B",
                  "ANPEP","APOE","CD14","CDKN1A","CHGA","CNTNAP2","EMP1",
                  "ERBB3","ERG","F5","FABP4","FASN","FGFR2","FOLH1","GDF15"
                  ,"GHR","GREM1","GSTP1","GUCY1A3","HPGD","HPN","IER3",
                  "IL1B","IL6","IL8","KRT15","KRT5","MME","MMP1","MMP7",
                  "MTHFD2","NAT1","NFKBIB","PDK4","PIK3R2","PLA2G7","PLAUR",
                  "PLN","PRAME","PSCA","RARRES1","SERPINB5","TIMP3","TRIM29"
                  )]
write.table(time0.28,"Time0.28.txt",col.names = T,row.names = F,sep="\t",quote = F)
# 49 genes

library(glasso)
t_static<- read.csv("t_static.csv",header=T)
t_static[,1]<-NULL
covA <- cov(t_static)
y <- seq(from=0.20, to=0.24, by=0.02) 
covR <- glassopath(covA,rholist=y) 
for (i in 1:3) {
  temp <- covR$wi[,,i] #wi:Estimated inverse covariance matrix, an array of dimension (nrow(s),ncol(n),length(rholist))
  write.table(x=temp, file=paste("t_static_matrix", i, sep = ""), sep="\t", row.names = F, col.names = F, quote = F)
}

# count non-zero number
library("Matrix")
a1<-read.table("t_static_matrix1",header=F) 
b1<-nnzero(a1,na.counted = TRUE)
print(b1-424) #0.2  670

a2<-read.table("t_static_matrix2",header=F) 
b2<-nnzero(a2,na.counted = TRUE)
print(b2-424) #0.22 414

a3<-read.table("t_static_matrix3",header=F) 
b3<-nnzero(a3,na.counted = TRUE)
print(b3-424) #0.24 264

colnames(a3)<-genes
rownames(a3)<-genes
genepair(a3,"t_static_0.24.txt")
uniqgene("t_static_0.24.txt","uniq_t_static_0.24.txt")

#abstract unique gene names from the files

a24<-read.table("t_static_0.24.txt",header=F)
a24_col2<-as.data.frame(a24[,2])
colnames(a24_col2)<-c("Gene ID")
b24<-unique(a24_col2)
b24t<-t(b24)
write.table(b24t,"b24t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE7701_10 GSE4196911_650.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,1:10]
t_time<-data.frame(t(time))
time0.24<-t_time[,c("ABCA5","ADIPOQ","ALOX12B","ALOX15","ALOX15B","ANPEP","ANXA2",
                    "APOE","BMP2","CAMKK2","CD14","CD24","CDKN1A","CHGA","CLASP2",
                    "CNTNAP2","CRP","CSF2","CYP3A4","EMP1","ERBB3","ERG","F5","FABP4",
                    "FASN","FGFR2","FOLH1","GDF15","GHR","GREB1","GREM1","GSTP1","GUCY1A3",
                    "HDAC9","HPGD","HPN","IER3","IHH","IL1B","IL6","IL8","INHBA","IRF7",
                    "ITGBL1","ITPR3","KDR","KHDRBS3","KRT15","KRT5","MEN1","MME","MMP1",
                    "MMP13","MMP7","MTHFD2","NAT1","NFATC4","NFKBIB","PDK4","PIK3R2","PLA2G7",
                    "PLAUR","PLN","PPFIA3","PRAME","PSCA","PTGDR","RARRES1","SERPINB5","SHH",
                    "SMO","SMOX","SULF1","TCL1B","TIMP3","TNFRSF13B","TRIM29","VDR","WIF1","ZNF10"
)]
write.table(time0.24,"Time0.24.txt",col.names = T,row.names = F,sep="\t",quote = F)
#80 genes

a34<-read.table("t_static_0.34.txt",header=F)
a34_col2<-as.data.frame(a34[,2])
colnames(a34_col2)<-c("Gene ID")
b34<-unique(a34_col2)
b34t<-t(b34)
write.table(b34t,"b34t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE7701_10 GSE4196911_650.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,1:10]
t_time<-data.frame(t(time))
time0.34<-t_time[,c("ABCA5","ALOX15","ALOX15B","ANPEP","CNTNAP2",
                    "ERG","F5","FASN","FGFR2","FOLH1","GDF15","GHR",
                    "HPGD","HPN","IL6","IL8","KRT15","KRT5","MMP7",
                    "NAT1","PIK3R2","PLA2G7","PLAUR","PLN","RARRES1",
                    "SERPINB5","TRIM29")]
write.table(time0.34,"Time0.34_27genes.txt",col.names = T,row.names = F,sep="\t",quote = F)
#27 genes