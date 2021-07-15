#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
library("ggplot2")
library("pheatmap")
library("ggfortify")
library("M3C")
library("Rtsne")
library("stringr")
library("reshape2")

setwd("/archive/data/hgc1074/yingchen/myiPS")

weighted_tpm <- readRDS("myips_tissue_1380.rds")
weighted_tpm[is.na(weighted_tpm)] <- 0
ips_tissue.anno <- as.data.frame(colnames(weighted_tpm))
colnames(ips_tissue.anno)<- "SampleID"
ips_tissue.anno$tissue<-str_split_fixed(ips_tissue.anno$SampleID, "_", 3)[,1]
ips_tissue.anno$tissue<-as.factor(ips_tissue.anno$tissue)
ips_tissue.anno$disease<-str_split_fixed(ips_tissue.anno$SampleID, "_", 3)[,2]
ips_tissue.anno$disease<-as.factor(ips_tissue.anno$disease)
ips_tissue.anno$disease <-gsub("IPS|iPS|KhES1\\.R2|KhES3\\.R2|KhES1\\.R1|KhES3\\.R1","healthy",ips_tissue.anno$disease) 
ips_tissue.anno$disease_detail <- str_split_fixed(ips_tissue.anno$SampleID, "_", 3)[,3]
ips_tissue.anno$disease_detail <- gsub("\\.R1|\\.R2","",ips_tissue.anno$disease_detail)
ips_tissue.anno$disease_detail <- gsub("1A|2A","healthy",ips_tissue.anno$disease_detail)

ips_tissue.anno$disease_detail <- gsub("^$","healthy",ips_tissue.anno$disease_detail)
ips_tissue.anno$tissue <- gsub(" ",".",ips_tissue.anno$tissue)
ips_tissue.anno$tissue <- as.factor(ips_tissue.anno$tissue)



disease_tissue <- as.data.frame(ips_tissue.anno$disease_detail[1:19])
disease_tissue_mat <- read.csv("disease_tissue_mat.csv",sep=",")
disease_tissue.ref <- cbind(disease_tissue,disease_tissue_mat)
disease_tissue.ref[is.na(disease_tissue.ref)]=0
rownames(disease_tissue.ref) <- disease_tissue.ref[,1]
#saveRDS(disease_tissue.ref,"disease_tissue.ref.rds")
colnames(disease_tissue.ref)[1]<-"disease_detail"
disease_tissue.long <- melt(disease_tissue.ref,id.vars="disease_detail")
colnames(disease_tissue.long) <- c("disease_detail","tissue","status")


ips_tissue.disease.anno <- merge(ips_tissue.anno,disease_tissue.long,by=c("disease_detail","tissue"),all.x=TRUE)
ips_tissue.disease.anno[is.na(ips_tissue.disease.anno)]<-0
rownames(ips_tissue.disease.anno) <- ips_tissue.disease.anno$SampleID
ips_tissue.disease.anno <- ips_tissue.disease.anno[,-3]

ips_tissue.disease.anno$status <- as.factor(ips_tissue.disease.anno$status)
#saveRDS(ips_tissue.disease.anno,"ips_tissue.disease.anno")

levels(ips_tissue.disease.anno$status) <-c("healthy","disease")

library("caret")
set.seed(999)
tpm_disease <- merge(t(weighted_tpm),ips_tissue.disease.anno, by="row.names",all.x=TRUE)
rownames(tpm_disease)<-tpm_disease$Row.names
tpm_disease <- tpm_disease[,-1]
tpm_disease <- tpm_disease[,-c(9334,9335,9336)]
intrain <- createDataPartition(y=tpm_disease$status,p=0.7,list=FALSE)
train_data <- tpm_disease[intrain,]
test_data <- tpm_disease[-intrain,]
anyNA(tpm_disease)
trainctrl <- trainControl(method="repeatedcv",number=10,repeats=3, classProbs=TRUE)

#svm_Linear <- train(status ~., data = train_data, method = "svmLinear",
#                 trControl=trainctrl,
#                 preProcess = c("center", "scale"),
#                 tuneLength = 10)
#svm_Linear
grid <- expand.grid(C=c(10^(-3),10^(-2),10^(-1),1,10,100,1000))
svm_Linear_Grid <- train(status ~.,data=train_data,method="svmLinear",trControl=trainctrl,tuneGrid=grid,preProcess=c("center","scale"),tunelength=10)


test_pred_grid <- predict(svm_Linear_Grid,newdata=test_data,type="prob")
save.image(file="1380ips.svm.prob.Rdata")
#ConfusionMatrix(test_pred_grid,test_data$status)

library("pROC")
rociPS <- roc(test_data$status,test_pred_grid$healthy)
pdf("ROC_test.pdf",paper="a4r")
ggroc(rociPS)
dev.off()



if(FALSE){
## PCA pic for 1380 samples
weighted_tpm <-as.data.frame(weighted_tpm)
weighted_tpm$var <- apply(weighted_tpm,1,var)
SelectGenes <- rownames(weighted_tpm[order(weighted_tpm$var,decreasing = T),][1:1000,])

pca_data <- t(weighted_tpm[,-ncol(weighted_tpm)][SelectGenes,])
pca_data <- log2(pca_data+1)
pca_data[is.na(pca_data)] <- 0

pcaResults <- prcomp(pca_data)

#pcaResults <- readRDS("gtex_pca_1000.rds")

pdf("ips.tissue.pca.pdf",paper="a4r")
autoplot(pcaResults,data=ips_tissue.anno,colour="tissue",label.show.legend = FALSE)
dev.off()
}

