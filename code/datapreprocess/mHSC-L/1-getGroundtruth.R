setwd("/home/bsb/mywork/rawdata/")
set.seed(123)
geneData <- read.csv("BEELINE-data/inputs/scRNA-Seq/mHSC-L/GeneOrdering.csv", row.names = 1)
exp=read.csv("BEELINE-data/inputs/scRNA-Seq/mHSC-L/ExpressionData.csv")
TF_list <- read.csv("mouse-tfs.csv")[,1]  
TF_list <- intersect(TF_list, rownames(geneData))
#Bonferroni 
geneData$P_adj <- p.adjust(geneData$VGAMpValue, method = "bonferroni")
sig <- subset(geneData,geneData$P_adj<0.01)
TF_sig <- intersect(rownames(sig), TF_list)
# get no-tf
nonTF_genes <- setdiff(rownames(sig), TF_sig)
nonTF_data <- sig[nonTF_genes,]
nonTF_data <- nonTF_data[order(nonTF_data$Variance,decreasing = T),]
#top500 top1000
top500=rownames(nonTF_data)[1:500]
max_rows <- min(1000, nrow(nonTF_data))  #632
top1000 <- rownames(nonTF_data)[1:max_rows]
combined500<- sample(c(TF_sig, top500))
combined1000<- sample(c(TF_sig, top1000))
#specific
speNet=read.csv("Networks/mouse/mHSC-ChIP-seq-network.csv")
speNet1000=subset(speNet,speNet$Gene1%in%combined1000 & speNet$Gene2%in%combined1000 &speNet$Gene1!=speNet$Gene2)
speNet500=subset(speNet,speNet$Gene1%in%combined500 & speNet$Gene2%in%combined500 &speNet$Gene1!=speNet$Gene2)
speTF1000 <- unique(speNet1000[,1])
speTarget1000 <- unique(speNet1000[,2])  
speTF500 <- unique(speNet500[,1])  
speTarget500 <- unique(speNet500[,2])  
#nonspecific 
nonspeNet=read.csv("Networks/mouse/Non-Specific-ChIP-seq-network.csv")
nonspeNet1000=subset(nonspeNet,nonspeNet$Gene1%in%combined1000 & nonspeNet$Gene2%in%combined1000  & nonspeNet$Gene1!=nonspeNet$Gene2)
nonspeNet500=subset(nonspeNet,nonspeNet$Gene1%in%combined500 & nonspeNet$Gene2%in%combined500& nonspeNet$Gene1!=nonspeNet$Gene2) 
nonspeTF1000 <- unique(nonspeNet1000[,1])  
nonspeTarget1000 <- unique(nonspeNet1000[,2])  
nonspeTF500 <- unique(nonspeNet500[,1]) 
nonspeTarget500 <- unique(nonspeNet500[,2])  
#string
stringNet=read.csv("Networks/mouse/STRING-network.csv")
stringNet1000=subset(stringNet,stringNet$Gene1%in%combined1000 & stringNet$Gene2%in%combined1000 & stringNet$Gene1!=stringNet$Gene2)
stringNet500=subset(stringNet,stringNet$Gene1%in%combined500 & stringNet$Gene2%in%combined500& stringNet$Gene1!=stringNet$Gene2) 
stringNet1000=unique(stringNet1000)
stringNet500=unique(stringNet500)
# stringNet1000 
stringTF1000 <- unique(stringNet1000[,1]) 
stringTarget1000 <- unique(stringNet1000[,2])  
stringTF500 <- unique(stringNet500[,1])  
stringTarget500 <- unique(stringNet500[,2])
#obtain exp
rownames(exp)=exp[,1]
exp=exp[,-1]
exp1000=exp[combined1000,]
exp500=exp[combined500,]
gene1000=data.frame(Gene=combined1000,index=0:(length(combined1000)-1))
gene500=data.frame(Gene=combined500,index=0:(length(combined500)-1))
rownames(gene500)=0:(length(combined500)-1)
rownames(gene1000)=0:(length(combined1000)-1)
tf1000=unique(c(TF_sig,speNet1000$Gene1,stringNet1000$Gene1,nonspeNet1000$Gene1))
tf1000=as.data.frame(tf1000)
tf1000 <- merge(tf1000, gene1000, by = 1)
tf1000 <- tf1000[order(tf1000$index), ]
tf500=unique(c(TF_sig,speNet500$Gene1,stringNet500$Gene1,nonspeNet500$Gene1))
tf500=as.data.frame(tf500)
tf500 <- merge(tf500, gene500, by = 1)
tf500 <- tf500[order(tf500$index), ]
rownames(tf1000)=0:(nrow(tf1000)-1)
colnames(tf1000)[1]="TF"
rownames(tf500)=0:(nrow(tf500)-1)
colnames(tf500)[1]="TF"
write.csv(exp1000,"../data/mHSC-L/exp1000.csv",quote=F,row.names = T)
write.csv(exp500,"../data/mHSC-L/exp500.csv",quote=F,row.names = T)
write.csv(gene1000,"../data/mHSC-L/Target1000.csv",quote = F)
write.csv(gene500,"../data/mHSC-L/Target500.csv",quote = F)
write.csv(tf1000,"../data/mHSC-L/TF1000.csv",quote = F)
write.csv(tf500,"../data/mHSC-L/TF500.csv",quote = F)
#spe
write.table(speTF1000, "../data/mHSC-L/SpecificDataset/groundtruth/TF1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(speTF500, "../data/mHSC-L/SpecificDataset/groundtruth/TF500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(speTarget1000, "../data/mHSC-L/SpecificDataset/groundtruth/Target1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(speTarget500, "../data/mHSC-L/SpecificDataset/groundtruth/Target500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(speNet1000, "../data/mHSC-L/SpecificDataset/groundtruth/net1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(speNet500, "../data/mHSC-L/SpecificDataset/groundtruth/net500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#nonspe
write.table(nonspeTF1000, "../data/mHSC-L/NonSpecificDataset/groundtruth/TF1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nonspeTF500, "../data/mHSC-L/NonSpecificDataset/groundtruth/TF500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nonspeTarget1000, "../data/mHSC-L/NonSpecificDataset/groundtruth/Target1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nonspeTarget500, "../data/mHSC-L/NonSpecificDataset/groundtruth/Target500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nonspeNet1000, "../data/mHSC-L/NonSpecificDataset/groundtruth/net1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(nonspeNet500, "../data/mHSC-L/NonSpecificDataset/groundtruth/net500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#string
write.table(stringTF1000, "../data/mHSC-L/StringDataset/groundtruth/TF1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(stringTF500, "../data/mHSC-L/StringDataset/groundtruth/TF500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(stringTarget1000, "../data/mHSC-L/StringDataset/groundtruth/Target1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(stringTarget500, "../data/mHSC-L/StringDataset/groundtruth/Target500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(stringNet1000, "../data/mHSC-L/StringDataset/groundtruth/net1000.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(stringNet500, "../data/mHSC-L/StringDataset/groundtruth/net500.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

