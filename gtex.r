#!/usr/local/package/r/3.6.0/bin/Rscript --vanilla
library("zFPKM")
library("ggplot2")
library("pheatmap")
library("ggfortify")
library("M3C")
library("Rtsne")
#setwd("/home/yingchen")

gtex_tpm <- readRDS("gtex_tpm.rds")
gtex_tpm <- as.data.frame(gtex_tpm)
gtex_tpm[is.na(gtex_tpm)] <- 0
gtex_tpm$sum <- rowSums(gtex_tpm)
gtex_cut <- subset(gtex_tpm,sum>0)
gtex_clean <- gtex_cut[,-ncol(gtex_cut)]
zTPM <- zFPKM(gtex_clean,assayName ="tpm")
activeGenes <- which(rowMeans(zTPM) > -3)
gtex_active <- gtex_clean[activeGenes,]
#saveRDS(gtex_active,"gtex_active.rds")


gtex_active <- readRDS("gtex_active.rds")
gtex_active$var <- apply(gtex_active,1,var)
gtex_active_clean<- gtex_active[,-ncol(gtex_active)]
#saveRDS(gtex_active,"gtex_var.rds")
selectGenes <- rownames(gtex_active[order(gtex_active$var,decreasing = T),][1:1000,])

gtex_anno <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",sep="\t",head=TRUE,quote = "",fill=NA)
gtex_anno$SAMPID <- gsub("-",".",gtex_anno$SAMPID)
gtex_anno <- subset(gtex_anno,gtex_anno$SAMPID %in% colnames(gtex_active))
gtex_anno_tissue <- cbind(gtex_anno$SAMPID, as.character(gtex_anno$SMTS))
colnames(gtex_anno_tissue) <- c("SampleID","Tissue")
gtex_anno_tissue <- as.data.frame(gtex_anno_tissue)
gtex_anno_tissue$Tissue <- as.factor(gtex_anno_tissue$Tissue)
rownames(gtex_anno_tissue) <- gtex_anno_tissue[,1]
#saveRDS(gtex_anno_tissue,"gtex_anno_tissue.rds")

#bar plot for tissue type. table for the tissue type

p_tissue <- ggplot(data=gtex_anno_tissue,aes(x=Tissue))+
	geom_bar(stat="count")+
	geom_text(stat="count",aes(label=..count..),size=2.5,vjust=-0.5)+
	theme_bw()+
	theme(axis.text=element_text(angle=45,vjust=0.5,hjust=1))
ggsave("tissue_counts_gg.pdf", p_tissue, width=8, height=6, units="in")
sample_info <- read.csv("gtex_sample_information.csv")



pca_data <- t(gtex_active_clean[selectGenes,])
pca_data <- log2(pca_data+1)
#pcaResults <- prcomp(pca_data)
#saveRDS(pcaResults,"gtex_pca_1000.rds")

#MDS plot from distance
MDSdist <- dist(pca_data)
#MDSresult <- cmdscale(MDSdistance,eig=TRUE, k=2)
pdf("gtex.tissue.MDS.pdf",paper="a4r")
autoplot(cmdscale(MDSdist, eig = TRUE),data=gtex_anno_tissue,colour="Tissue",label.show.legend=FALSE)
dev.off()
if(FALSE){
#pcaResults <- readRDS("gtex_pca_1000.rds")

pdf("gtex.tissue.pca.pdf",paper="a4r")
autoplot(pcaResults,data=gtex_anno_tissue,colour="Tissue",label.show.legend = FALSE)
dev.off()


#tsne
set.seed(999)
gtex_active <- readRDS("gtex_active.rds")
gtex_anno_tissue <- readRDS("gtex_anno_tissue.rds")
colnames(gtex_anno_tissue) <- c("ID","Tissue")
gtex_active$var <- apply(gtex_active,1,var)
gtex_active_clean<- gtex_active[,-ncol(gtex_active)]
selectGenes <- rownames(gtex_active[order(gtex_active$var,decreasing = T),][1:1000,])
tsne_data <- gtex_active_clean[selectGenes,]
tsne_data <- log2(tsne_data+1)
saveRDS(tsne_data,"tsne_data.rds")

#res <- M3C(tsne_data, des=gtex_anno_tissue, removeplots=TRUE, iters=25,fsize=8,maxK=20,analysistype="chi")
#saveRDS(res,"M3Cresults.rds")

#res <- M3C(tsne_data, des=gtex_anno_tissue,method=2,iters=25,maxK=20,analysistype="chi")
#saveRDS(res,"M3Cresults_fast.rds")
tsne <- Rtsne(t(tsne_data),dims=2,perplexity=30,verbose=TRUE,max_iter=500,num_threads=0)
saveRDS(tsne,"rtsne_result.rds")

embedding <- as.data.frame(tsne$Y)
embedding$Tissue <- as.factor(gtex_anno_tissue$Tissue)

p <- ggplot(embedding, aes(x=V1, y=V2, color=Tissue)) +
     geom_point(size=1.25) +
     guides(colour = guide_legend(override.aes = list(size=6))) +
     xlab("") + ylab("") +
     ggtitle("t-SNE Gtex_only") +
     theme_light(base_size=20) +
     theme(strip.background = element_blank(),
           strip.text.x     = element_blank(),
           axis.text.x      = element_blank(),
           axis.text.y      = element_blank(),
           axis.ticks       = element_blank(),
           axis.line        = element_blank(),
           panel.border     = element_blank())

ggsave("rtsne_result_gg.pdf", p, width=8, height=6, units="in")


gender <- c("Female","Male")
gender_counts <-c(11584,5798)

table_gender <- data.frame(gender,gender_counts)
colnames(table_gender) <-c("Gender","Counts")
table_gender$prop <- table_gender$Counts / 17382
table_gender$lab.ypos = cumsum(table_gender$prop) - 0.5*table_gender$prop

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

p_gender <- ggplot(table_gender, aes(x = "", y = prop, fill = Gender)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()
ggsave("gtex_gender.pdf", p_gender,width=8, height=6,units="in")
# Gender Counts      prop  lab.ypos
#1 Female  11584 0.6664365 0.3332183
#2   Male   5798 0.3335635 0.8332183
age <- c("20-29","30-39","40-49","50-59","60-69","70-79")
age_counts <-c(1320,1323,2702,5615,5821,601)
table_age <-data.frame(age,age_counts)

colnames(table_age) <-c("Age","Counts")
table_age$prop <- table_age$Counts / 17382
table_age$lab.ypos = cumsum(table_age$prop) - 0.5*table_age$prop

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

p_age <- ggplot(table_age, aes(x = "", y = prop, fill = Age)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void()
ggsave("gtex_age.pdf", p_age,width=8, height=6,units="in")




#colors <- rainbow(length(unique(gtex_anno_tissue$Tissue)))
#names(colors) <- unique(gtex_anno_tissue$Tissue)
#pdf("gtex_rtsne.pdf",paper="a4r")
#par(mgp=c(2.5,1,0))
#plot(tsne$Y, t='n', main="tSNE",xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2,"cex.lab"=1.5)
#text(tsne$Y, labels=gtex_anno_tissue$Tissue,col=colors[gtex_anno_tissue$Tissue])
#dev.off()


#pdf("gtex.tissue.tsne.pdf",paper="a4r")
#tsne(tsne_data)
#dev.off()

}
