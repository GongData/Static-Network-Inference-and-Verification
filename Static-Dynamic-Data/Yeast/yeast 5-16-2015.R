library(Biobase)
library(GEOquery)
library(limma)
gset<-getGEO("GSE19213", GSEMatrix =TRUE)
gset<-gset[[1]]
show(gset)
str(gset)
e<-data.frame(exprs(gset))
boxplot(log(e))

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# get top 2000 different expressed genes between control and H2O2 treatment
# group names for all samples
sml <- c("G0","G0","G0","G0","X","X","X","X","G1","G1","G1","G1","X","X","X","X");

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=5000)
tT<-tT[-1,]
tTID <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))

  
  
#gene symbol
gene<-as.character(gset@featureData$Gene.Symbol)
ID<-as.character(gset@featureData$ID)
IDs<-data.frame(cbind(ID,gene))

rows<-rownames(tTID)
# loop over the array to create a vector of gene symbols
length(rows)
symb<-rep(0,length(rows))
for (i in 1:length(rows)){
  symb[i]<-as.character(IDs[which(IDs[,1]==rows[i]),2])
}

# We don't actually have gene symbols for all rows, we need to address for it
sel<-apply(as.matrix(symb,ncol=1),1,
           function(x){if(x=="")return(FALSE)else return(TRUE)})

selexpr<-c()
selsymb<-c()

for (i in 1:nrow(tTID)){
  if (sel[i]){
    selexpr<-rbind(selexpr,tTID[i,])
    selsymb<-rbind(selsymb,symb[i])
  }
}
selexpr[,"Gene Symbol"]<-selsymb
write.table(x =selexpr, file = 'GSE19213Top.csv', sep =',',col.names=T, row.names=F, quote = F, na = '')

#get the expression value of the top 2000 diferent expressed genes
H2O2<-e
y<-c("GSM476089","GSM476090","GSM476091","GSM476092")
H2O2<-H2O2[y]
Top<-read.csv('GSE19213Top.csv',header=T)
m1 <- merge(Top, H2O2, by.x = "ID", by.y = "row.names")

static<-m1
static[,2]<-NULL #repeat 5 times
symb<-as.character(static[,2])
static[,2]<-NULL
rows<-static[,1]
static[,1]<-NULL
rownames(static)<-rows

# get average value for each gene
ave<-aggregate(static, by=list(symb), 
               FUN=mean)
ave<-ave[-1,]
rows<-ave[,1]
ave[,1]<-NULL
rownames(ave)<-rows
time<-read.csv("ave_GSE62120H2O2_time.csv",header=T)
rows<-time[,1]
time[,1]<-NULL
rownames(time)<-rows

#conbime time series and static data together
m<-merge(ave,time,by.x="row.names",by.y="row.names")
write.table(x =m, file = 'GSE19213_4 GSE6212_11.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')

Total<-read.csv('GSE19213_4 GSE6212_11.csv',header=T)
rows<-Total[,1]
Total[,1]<-NULL
rownames(Total)<-rows

static<-Total[1:4]
t_static<-data.frame(t(static))
write.table(x =t_static, file = 't_static.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')

library(glasso)
t_static<- read.csv("t_static.csv",header=T)
t_static[,1]<-NULL
covA <- cov(t_static)
y <- seq(from=0.1, to=0.2, by=0.02) 
covR <- glassopath(covA,rholist=y) 
for (i in 1:6) {
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


genes<-rownames(static)
colnames(a1)<-genes
rownames(a1)<-genes

colnames(a2)<-genes
rownames(a2)<-genes

colnames(a3)<-genes
rownames(a3)<-genes

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

genepair(a1,"t_static_0.26.txt")
uniqgene("t_static_0.26.txt","uniq_t_static_0.26.txt")

genepair(a2,"t_static_0.30.txt")
uniqgene("t_static_0.30.txt","uniq_t_static_0.30.txt")

genepair(a3,"t_static_0.38.txt")
uniqgene("t_static_0.38.txt","uniq_t_static_0.38.txt")

#abstract unique gene names from the files

a26<-read.table("t_static_0.26.txt",header=F)
a26_col2<-as.data.frame(a26[,2])
colnames(a26_col2)<-c("Gene ID")
b26<-unique(a26_col2)
b26t<-t(b26)
write.table(b26t,"b26t.txt",col.names = F,row.names = T,sep=",",quote = T)

m<-read.csv('GSE19213_4 GSE6212_11.csv',header=T)
rows<-m[,1]
m[,1]<-NULL
rownames(m)<-rows
time<-m[,5:15]
t_time<-data.frame(t(time))
time0.26<-t_time[,c("ABP1","ADA2","ADE4","AIR1","ARX1","ATG17",
                    "CGR1","CNS1","COG3","CYR1","DAD2","DAK2",
                    "DBP3","DBP7","DBP8","DIM1","DRS1","EBP2",
                    "EPL1","FBP1","FMT1","GIT1","HOP1","HOS2",
                    "KAP123","MED7","MRPL4","MSA1","MSH1","MSH3",
                    "MTR4","NRM1","PAN1","PDS5","PHO4","POL4","PPR1",
                    "PUF3","PUF6","PUS2","PUT4","REC8","RNH201","RPA12",
                    "RPA43","RRP12","SCM3","SDS23","SLD5","SNU66","SPB1",
                    "SPO20","SWF1","TEA1","TRM9","UTP10","UTP13","UTP20",
                    "UTP4","VRP1","YOX1","ZDS1"

)]
write.table(time0.26,"Time0.26.txt",col.names = T,row.names = F,sep="\t",quote = F)
# 62 genes

a30<-read.table("t_static_0.30.txt",header=F)
a30_col2<-as.data.frame(a30[,2])
colnames(a30_col2)<-c("Gene ID")
b30<-unique(a30_col2)
b30t<-t(b30)
write.table(b30t,"b30t.txt",col.names = F,row.names = T,sep=",",quote = T)

time0.30<-t_time[,c("ADA2","AIR1","ARX1","ATG17","CGR1","CNS1",
                    "DAD2","DAK2","DBP7","DIM1","DRS1","FMT1",
                    "GIT1","HOP1","HOS2","MRPL4","MSA1","MSH1",
                    "MSH3","MTR4","PAN1","PHO4","POL4","PPR1",
                    "PUT4","RPA12","RPA43","RRP12","SNU66","SPB1",
                    "SPO20","SWF1","TEA1","UTP13","VRP1","YOX1","ZDS1"                    
)]
write.table(time0.30,"Time0.30.txt",col.names = T,row.names = F,sep="\t",quote = F)
#37 genes

a38<-read.table("t_static_0.38.txt",header=F)
a38_col2<-as.data.frame(a38[,2])
colnames(a38_col2)<-c("Gene ID")
b38<-unique(a38_col2)
b38t<-t(b38)
write.table(b38t,"b38t.txt",col.names = F,row.names = T,sep=",",quote = T)

time0.38<-t_time[,c("ADA2","AIR1","ARX1","CGR1","DAK2","FMT1","GIT1",
                    "HOS2","MRPL4","MTR4","PHO4","POL4","RPA12","RRP12",
                    "SNU66","SPB1","SWF1","TEA1","VRP1","YOX1","ZDS1"
                  )]
write.table(time0.38,"Time0.38.txt",col.names = T,row.names = F,sep="\t",quote = F)
#21 genes



