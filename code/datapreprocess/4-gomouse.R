library(HGNChelper)
library(clusterProfiler)
library(org.Mm.eg.db)
library(proxy)
set.seed(123)
calculate_go_similarity <- function(cell_type, net_type, gene_type, gene_size) {
  input_path <- file.path("/home/bsb/mywork/data", cell_type, net_type, "groundtruth")
  gene_file <- paste0(gene_type, gene_size, ".txt")
  gene_list <- read.table(file.path(input_path, gene_file))[, 1]
  final_mat <- matrix(NA, 
                      nrow = length(gene_list), 
                      ncol = length(gene_list),
                      dimnames = list(gene_list, gene_list))
  corrected <- checkGeneSymbols(gene_list, species = "mouse")
  go_mapping <- bitr(corrected$Suggested.Symbol, "SYMBOL", "GO", org.Mm.eg.db)
  # 
  merged <- merge(corrected, go_mapping, by.x=3, by.y=1)
  go_mat <- +(table(merged[,2], merged[,4]) > 0)
  anno_sim <- as.matrix(simil(go_mat, method = "cosine"))
  diag(anno_sim) <- 1
  final_mat[rownames(anno_sim), rownames(anno_sim)] <- anno_sim
  #missing fill
  na_mask <- is.na(final_mat)
  valid_vals <- final_mat[!na_mask & final_mat != 1]
  final_mat[na_mask] <- sample(valid_vals, size = sum(na_mask), replace = TRUE)
  diag(final_mat) <- 1
  colnames(final_mat)[1] <- paste("\t", colnames(final_mat)[1], sep = "")
  prefix <- if (gene_type == "TF") {
    "tf"
  } else if (gene_type == "Target") {
    "tg"
  } else {
    gene_type  
  }
  output_dir <- file.path("/home/bsb/mywork/data", cell_type, net_type, "similarity")
  output_file <- file.path(output_dir, paste0(prefix, gene_size, "gocos.txt"))
  write.table(final_mat, output_file, quote = F, row.names = T, col.names = T, sep = "\t")
  missing_genes <- setdiff(gene_list, rownames(anno_sim))  
  cat(sprintf("celltype:%s | net:%s | %s%s | gene num:%d | missing num:%d\n",
              cell_type, net_type, gene_type, gene_size, length(gene_list), length(missing_genes)))
}

cat("SpecificDataset:\n")
calculate_go_similarity("mESC", "SpecificDataset", "TF", 500)
calculate_go_similarity("mESC", "SpecificDataset", "TF", 1000)
calculate_go_similarity("mESC", "SpecificDataset", "Target", 500)
calculate_go_similarity("mESC", "SpecificDataset", "Target", 1000)
cat("\n NonSpecificDataset:\n")
calculate_go_similarity("mESC", "NonSpecificDataset", "TF", 500)
calculate_go_similarity("mESC", "NonSpecificDataset", "TF", 1000)
calculate_go_similarity("mESC", "NonSpecificDataset", "Target", 500)
calculate_go_similarity("mESC", "NonSpecificDataset", "Target", 1000)
cat("\n StringDataset:\n")
calculate_go_similarity("mESC", "StringDataset", "TF", 500)
calculate_go_similarity("mESC", "StringDataset", "TF", 1000)
calculate_go_similarity("mESC", "StringDataset", "Target", 500)
calculate_go_similarity("mESC", "StringDataset", "Target", 1000)

cat("SpecificDataset:\n")
calculate_go_similarity("mDC", "SpecificDataset", "TF", 500)
calculate_go_similarity("mDC", "SpecificDataset", "TF", 1000)
calculate_go_similarity("mDC", "SpecificDataset", "Target", 500)
calculate_go_similarity("mDC", "SpecificDataset", "Target", 1000)
cat("\n NonSpecificDataset:\n")
calculate_go_similarity("mDC", "NonSpecificDataset", "TF", 500)
calculate_go_similarity("mDC", "NonSpecificDataset", "TF", 1000)
calculate_go_similarity("mDC", "NonSpecificDataset", "Target", 500)
calculate_go_similarity("mDC", "NonSpecificDataset", "Target", 1000)
cat("\n StringDataset:\n")
calculate_go_similarity("mDC", "StringDataset", "TF", 500)
calculate_go_similarity("mDC", "StringDataset", "TF", 1000)
calculate_go_similarity("mDC", "StringDataset", "Target", 500)
calculate_go_similarity("mDC", "StringDataset", "Target", 1000)

cat("SpecificDataset:\n")
calculate_go_similarity("mHSC-E", "SpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-E", "SpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-E", "SpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-E", "SpecificDataset", "Target", 1000)
cat("\n NonSpecificDataset:\n")
calculate_go_similarity("mHSC-E", "NonSpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-E", "NonSpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-E", "NonSpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-E", "NonSpecificDataset", "Target", 1000)
cat("\n StringDataset:\n")
calculate_go_similarity("mHSC-E", "StringDataset", "TF", 500)
calculate_go_similarity("mHSC-E", "StringDataset", "TF", 1000)
calculate_go_similarity("mHSC-E", "StringDataset", "Target", 500)
calculate_go_similarity("mHSC-E", "StringDataset", "Target", 1000)

cat("SpecificDataset:\n")
calculate_go_similarity("mHSC-GM", "SpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-GM", "SpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-GM", "SpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-GM", "SpecificDataset", "Target", 1000)
cat("\n NonSpecificDataset:\n")
calculate_go_similarity("mHSC-GM", "NonSpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-GM", "NonSpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-GM", "NonSpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-GM", "NonSpecificDataset", "Target", 1000)
cat("\n StringDataset:\n")
calculate_go_similarity("mHSC-GM", "StringDataset", "TF", 500)
calculate_go_similarity("mHSC-GM", "StringDataset", "TF", 1000)
calculate_go_similarity("mHSC-GM", "StringDataset", "Target", 500)
calculate_go_similarity("mHSC-GM", "StringDataset", "Target", 1000)

cat("SpecificDataset:\n")
calculate_go_similarity("mHSC-L", "SpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-L", "SpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-L", "SpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-L", "SpecificDataset", "Target", 1000)
cat("\n NonSpecificDataset:\n")
calculate_go_similarity("mHSC-L", "NonSpecificDataset", "TF", 500)
calculate_go_similarity("mHSC-L", "NonSpecificDataset", "TF", 1000)
calculate_go_similarity("mHSC-L", "NonSpecificDataset", "Target", 500)
calculate_go_similarity("mHSC-L", "NonSpecificDataset", "Target", 1000)
cat("\n StringDataset:\n")
calculate_go_similarity("mHSC-L", "StringDataset", "TF", 500)
calculate_go_similarity("mHSC-L", "StringDataset", "TF", 1000)
calculate_go_similarity("mHSC-L", "StringDataset", "Target", 500)
calculate_go_similarity("mHSC-L", "StringDataset", "Target", 1000)


