
##################################################################
#############################s############### COAD ############################
##################################################################


Raw.Data.Path<-"E:/Dropbox/Work_temp/IO_nanoStrings_CRC/TCGA"
setwd(Raw.Data.Path)

COAD.Clinical<-read.table(file="survival_COAD_survival.txt",sep="\t",head=TRUE,quote = "\"")
head(COAD.Clinical)

COAD.Clinical2<-read.table(file="TCGA.COAD.sampleMap_COAD_clinicalMatrix",sep="\t",head=TRUE,quote = "\"")
head(COAD.Clinical2[,1:10])

TCGA.temp2<-read.table(file="TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp2[,1:10])
TCGA.temp<-read.table(file="TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp[,1:10])
dim(TCGA.temp)
dim(TCGA.temp2)

# TCGA.temp2[,-1]<-head(TCGA.temp2[,1:10])
# match(colnames(TCGA.temp),"TCGA.AG.3611.01")

# matrix.temp<-2^TCGA.temp[,-1]-1
# TCGA.temp2<-log2(matrix.temp*10^6/t(matrix(colSums(matrix.temp),512,60483))+1)
# head(TCGA.temp2[,1:10])

COAD.T.index<-match(unique(paste0(COAD.Clinical$X_PATIENT,"-01")),COAD.Clinical$sample,nomatch=0)
sum(COAD.T.index!=0)
COAD.N.index<-match(unique(paste0(COAD.Clinical$X_PATIENT,"-11")),COAD.Clinical$sample,nomatch=0)
sum(COAD.N.index!=0)

# COAD.T.index<-match(paste0(colnames(TCGA.temp2),"-01"),COAD.Clinical$sample,nomatch=0)
# sum(COAD.T.index!=0)

COAD.matchT.index<-match(unique(paste0(COAD.Clinical$X_PATIENT[COAD.N.index],"-01")),COAD.Clinical$sample,nomatch=0)
cbind(COAD.Clinical$sample[COAD.matchT.index],COAD.Clinical$sample[COAD.N.index])

COAD.Clinical.T<-COAD.Clinical[COAD.T.index,]

common.names<-intersect(substr(colnames(TCGA.temp),1,15),make.names(COAD.Clinical.T$sample))
common.names
length(common.names)

common.names2<-intersect(substr(colnames(TCGA.temp2),1,15),make.names(COAD.Clinical.T$sample))
common.names2
length(common.names2)

common.namesN<-intersect(substr(colnames(TCGA.temp2),1,15),make.names(COAD.Clinical$sample[COAD.N.index]))
common.namesN
length(common.namesN)


