library(dplyr)
library(data.table) 
options(stringsAsFactors = F)
#1.1 load data
human_struc = read.table("/home/bsb/mywork/data/middlefile/humanresult/tf_similarity.tsv", 
                         header=F, sep="\t", stringsAsFactors=F)
mouse_struc = read.table("/home/bsb/mywork/data/middlefile/mouseresult/tf_similarity.tsv", 
                         header=F, sep="\t", stringsAsFactors=F)
colnames(human_struc) = c("query", "target", "alntmscore", "prob", "LDDT","qtmscore","ttmscore","evalue")
colnames(mouse_struc)=colnames(human_struc)
#1.2 get uniprot id
extract_uniprot_id <- function(af_id) {
  strsplit(af_id, "-")[[1]][2]
}
human_struc$query <- sapply(human_struc$query, extract_uniprot_id)
human_struc$target <- sapply(human_struc$target, extract_uniprot_id)
mouse_struc$query <- sapply(mouse_struc$query, extract_uniprot_id)
mouse_struc$target <- sapply(mouse_struc$target, extract_uniprot_id)
mouse_struc <- mouse_struc %>% 
  group_by(query, target) %>% 
  summarise(
    alntmscore = mean(alntmscore, na.rm = TRUE),
    prob = mean(prob, na.rm = TRUE),
    LDDT = mean(LDDT, na.rm = TRUE),
    qtmscore = mean(qtmscore, na.rm = TRUE),
    ttmscore = mean(ttmscore, na.rm = TRUE),
    evalue = mean(evalue, na.rm = TRUE),
    .groups = 'drop'  
  )
#1.3 id change
humanid=read.table("/home/bsb/mywork/data/middlefile/human_tf_uniprot.txt",header = T)
mouseid=read.table("/home/bsb/mywork/data/middlefile/mouse_tf_uniprot.txt",header = T)
human_struc=merge(human_struc,humanid,by.x=1,by.y=2)[,c(9,2:8)]
colnames(human_struc)[1]="query"
human_struc=merge(human_struc,humanid,by.x=2,by.y=2)[,c(2,9,3:8)]
colnames(human_struc)[2]="target"
mouse_struc=merge(mouse_struc,mouseid,by.x=1,by.y=2)[,c(9,2:8)]
colnames(mouse_struc)[1]="query"
mouse_struc=merge(mouse_struc,mouseid,by.x=2,by.y=2)[,c(2,9,3:8)]
colnames(mouse_struc)[2]="target"
#1.4 add mean max
human_struc$mean <- (human_struc$qtmscore + human_struc$ttmscore)/2
human_struc$max  <- pmax(human_struc$qtmscore, human_struc$ttmscore, na.rm = TRUE)
mouse_struc$mean <- (mouse_struc$qtmscore + mouse_struc$ttmscore)/2
mouse_struc$max  <- pmax(mouse_struc$qtmscore, mouse_struc$ttmscore, na.rm = TRUE)

#2 param
human_celltype_list <- c("hESC", "hHep")  
mouse_celltype_list <- c("mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")
net_list        <- c("SpecificDataset", "NonSpecificDataset", "StringDataset") 
size_list     <- c(500, 1000) 
score_cols       <- c("max")# "alntmscore","prob","LDDT","qtmscore","ttmscore","mean",

#3.1 human full matrix
cat("Start computing full human TF similarity matrix")
setDT(human_struc) 
human_all_matrix <- list()
for (score_col in score_cols) {
  mat <- dcast(human_struc, query ~ target, value.var = score_col, fill = NA)
  # Convert to matrix format, row names = query, column names = target
  mat_matrix <- as.matrix(mat[, -1]) 
  rownames(mat_matrix) <- mat$query
  human_all_matrix[[score_col]] <- mat_matrix
}
#  3.2 Mouse full matrix
cat("Start precomputing full mouse TF similarity matrix ")
setDT(mouse_struc) 
mouse_all_matrix <- list()
for (score_col in score_cols) {
  mat <- dcast(mouse_struc, query ~ target, value.var = score_col, fill = NA)
  mat_matrix <- as.matrix(mat[, -1])
  rownames(mat_matrix) <- mat$query
  mouse_all_matrix[[score_col]] <- mat_matrix
}
#4.1 human
generate_tf_similarity_matrix <- function(celltype, net_type, tf_size, all_matrix) {
  tf_file_path  <- paste0("/home/bsb/mywork/data/", celltype, "/", net_type, "/groundtruth/TF",tf_size,".txt")
  save_dir_path <- paste0("/home/bsb/mywork/data/", celltype, "/", net_type, "/similarity/")
  tf_gene_list  <- read.table(tf_file_path, header = F)[,1]
  for (score_col in score_cols) {
    score_matrix <- all_matrix[[score_col]][tf_gene_list, tf_gene_list]
    save_file_name <- paste0(save_dir_path, "tf", tf_size, "_", score_col, ".txt")
    colnames(score_matrix)[1] <- paste0("\t", colnames(score_matrix)[1])
    write.table(score_matrix, file = save_file_name, quote = F, row.names = T, col.names = T, sep = "\t")
  }
}
#
cat("===== hESC/hHep TF similarity analysis start =====\n")
for(celltype in human_celltype_list){
  for(net in net_list){
    for(size in size_list){
      generate_tf_similarity_matrix(celltype, net, size, human_all_matrix)
    }
  }
}
cat("===== mDC/mESC/mHSC series TF similarity analysis start =====\n")
for(celltype in mouse_celltype_list){
  for(net in net_list){
    for(size in size_list){
      generate_tf_similarity_matrix(celltype, net, size, mouse_all_matrix)
    }
  }
}
cat("\n finish all")

