library(Biobase)
library(GEOquery)
library(limma)
gset<-getGEO("GSE62120", GSEMatrix =TRUE)
gset<-gset[[1]]
show(gset)
str(gset)
e<-data.frame(exprs(gset))
boxplot(log(e))
genes<-gset@featureData$"Gene Symbol"
length(genes)

platf<-getGEO(annotation(gset),AnnotGPL=TRUE)
str(platf)
IDs<-attr(dataTable(platf),"table")[,c("ID","Gene symbol")]
head(IDs)

# creat a new column called Gene_SYMBOL which indicate the gene symbol of 
# the probles ID in the data.frame
rows<-rownames(e)
head(rows)
# loop over the array to create a vector of gene symbols
length(rows)
symb<-rep(0,length(rows))
for (i in 1:length(rows)){
  symb[i]<-as.character(IDs[which(IDs[,1]==rows[i]),2])
}
symb[1:10]

# get average value for each gene
ave<-aggregate(e, by=list(symb), 
               FUN=mean)
ave<-ave[-1,]
rows<-ave[,1]
ave[,1]<-NULL
rownames(ave)<-rows
y<-c("GSM1519986","GSM1519987","GSM1519988","GSM1519989","GSM1519990","GSM1519991","GSM1519992",
     "GSM1519993","GSM1519994","GSM1519995","GSM1519996")
time<-ave[y]
write.table(x =time, file = 'ave_GSE62120H2O2_time.csv', sep =',',col.names=NA, row.names=T, quote = F, na = '')



