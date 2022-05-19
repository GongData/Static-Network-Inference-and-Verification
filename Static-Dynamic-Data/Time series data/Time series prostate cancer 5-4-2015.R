# study design of the source data:
# Gene Expression Profiling of Prostate Benign and Malignant Tissue
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55945

library(Biobase)
library(GEOquery)
# static data
gset55945<-getGEO("GSE55945", GSEMatrix =TRUE)
gset55945<-gset55945[[1]]
show(gset55945)
str(gset55945)


#check if the dataset appears to be normalized
e2<-exprs(gset55945)
show(pData(phenoData(gset55945))[,1:3])
boxplot(log(e2))
# Already normolized but did not get log transformation

# time series data
gset30174<-getGEO("GSE30174", GSEMatrix =TRUE)
gset30174<-gset30174[[1]]
show(gset30174)
str(gset30174)


#check if the dataset appears to be normalized
e<-exprs(gset30174)
show(pData(phenoData(gset30174))[,1:3])
boxplot(log(e))
# Did not normolized

# work on the gset55945 firstly
str(gset55945)
symb<-gset55945@featureData @data $"Gene Symbol"  
symb<-as.character(symb)
# We don't actually have gene symbols for all rows, we need to address for it
sel<-apply(as.matrix(symb,ncol=1),1,
           function(x){if(x=="")return(FALSE)else return(TRUE)})
levels(as.factor(sel))
selexpr<-c()
selsymb<-c()
for (i in 1:nrow(e2)){
  if (sel[i]){
    selexpr<-rbind(selexpr,e2[i,])
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


# check the expression values
min(ave)
max(ave) #19160.83 too high
summary(ave)
# delete genes with expression value more thant 10000.00
ave_sub<-subset(ave, GSM1348933<10000 &GSM1348934<10000 &GSM1348935<10000
                &GSM1348936<10000 &GSM1348937<10000 &GSM1348938<10000 
                &GSM1348939<10000 &GSM1348940<10000 &GSM1348941<10000 
                &GSM1348942<10000 &GSM1348943<10000 &GSM1348944<10000
                &GSM1348945<10000 &GSM1348946<10000 &GSM1348947<10000
                &GSM1348948<10000 &GSM1348949<10000 &GSM1348950<10000
                &GSM1348951<10000 &GSM1348952<10000 &GSM1348953<10000)
summary(ave_sub)

# log2 transform
qx <- as.numeric(quantile(ave_sub, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ave_sub[which(ave_sub <= 0)] <- NaN
            ave_sub_l <- log2(ave_sub) }
summary(ave_sub_l)
a<-read.csv("ave_GSE55945.csv",header=T)
rows<-a[,1]
a[,1]<-NULL
rownames(a)<-rows
source_name<-a[21050,]
source_name<-as.factor(source_name)
ave_sub_l<-rbind(ave_sub_l,source_name)
ave_sub_l_2<-ave_sub_l[1:20989,]


# compare the malignant tissue /  benign tissue
library(limma)
design<- model.matrix(~source_name+0,t_ave_sub_l)
colnames(design) <- c("benign","malignant")
fit <- lmFit(ave_sub_l_2, design)




#subgroup the data: malignant tissue /  benign tissue
summary(gset55945$source_name_ch1)
gset5_mal<-gset55945[,which(gset55945$source_name_ch1=="malignant tissue_TMPRSS2:ERG fusion negative"
                            |gset55945$source_name_ch1=="malignant tissue_TMPRSS2:ERG fusion positive")]
summary(gset5_mal$source_name_ch1)
gset5_beni<-gset55945[,which(gset55945$source_name_ch1=="benign tissue")]
summary(gset5_beni$source_name_ch1)

# We don't actually have gene symbols for all rows, we need to address for it
symb<-gset55945@featureData @data $"Gene Symbol"  
symb<-as.character(symb)
sel<-apply(as.matrix(symb,ncol=1),1,
           function(x){if(x=="")return(FALSE)else return(TRUE)})
levels(as.factor(sel))
selexpr_mal<-c()
selsymb_mal<-c()
selexpr_beni<-c()
selsymb_beni<-c()
e_mal<-exprs(gset5_mal)
e_beni<-exprs(gset5_beni)
for (i in 1:nrow(e_mal)){
  if (sel[i]){
    selexpr_mal<-rbind(selexpr_mal,e_mal[i,])
    selsymb_mal<-rbind(selsymb_mal,symb[i])
  }
}
selexpr_mal<-data.frame(selexpr_mal)
for (i in 1:nrow(e_beni)){
  if (sel[i]){
    selexpr_beni<-rbind(selexpr_beni,e_beni[i,])
    selsymb_beni<-rbind(selsymb_beni,symb[i])
  }
}
selexpr_beni<-data.frame(selexpr_beni)

# get average value for each gene

ave_mal<-aggregate(selexpr_mal, by=list(selsymb_mal), 
               FUN=mean)
write.table(x =ave_mal, file ="ave_GSE55945_mal.csv", sep =',', col.names=T,row.names=F, quote = F, na = '')
ave_beni<-aggregate(selexpr_beni, by=list(selsymb_beni), 
                   FUN=mean)
write.table(x =ave_beni, file ="ave_GSE55945_beni.csv", sep =',', col.names=T,row.names=F, quote = F, na = '')




