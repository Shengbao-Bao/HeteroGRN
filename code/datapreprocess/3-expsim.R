set.seed(1)
cell_types <- c("hESC","hHep","mHSC-GM", "mHSC-E","mHSC-L", "mESC", "mDC")
base_path <- "/home/bsb/mywork/data/"
sample_sizes <- c(500, 1000)
gene_types <- c("tf", "tg")
dataset_types <- c("NonSpecific", "Specific", "String")
calc_pearson_similarity <- function(cell_type, gene_type, sample_size, dataset_type, base_path) {
  main_dir <- paste0(base_path, cell_type, "/", dataset_type, "Dataset/")
  gene_file <- paste0(main_dir, "groundtruth/", 
                      ifelse(gene_type=="tf", "TF", "Target"), sample_size, ".txt")
  exp_file <- paste0(main_dir, "../exp", sample_size, ".csv")  
  save_dir <- paste0(main_dir, "similarity/")
  #load data
  genes <- scan(gene_file, what = character())
  exp_data <- read.csv(exp_file, row.names = 1)
  gene_exp <- exp_data[genes,]  
  message("process：", cell_type, "-", dataset_type, "-", gene_type, sample_size, 
          " | valid gene num：", nrow(gene_exp))
  exp_t <- t(gene_exp)  
  pearson <- abs(cor(exp_t, method = "pearson"))  
  cosine <- lsa::cosine(exp_t)
  
  fix_colname <- function(mat) {
    if (ncol(mat) > 0) colnames(mat)[1] <- paste0("\t", colnames(mat)[1])
    mat
  }
  # 
  write.table(fix_colname(gene_exp), 
              paste0(save_dir, gene_type, sample_size, "exp.txt"), 
              row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(fix_colname(pearson), 
              paste0(save_dir, gene_type, sample_size, "pearson.txt"), 
              row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  message("save：", paste0(save_dir, gene_type, sample_size, "pearson.txt"), "\n")
}

for (cell in cell_types) {
  for (dataset in dataset_types) {
    for (gene in gene_types) {
      for (size in sample_sizes) {
        tryCatch({
          calc_pearson_similarity(
            cell_type = cell,
            gene_type = gene,
            sample_size = size,
            dataset_type = dataset,
            base_path = base_path
          )
        }, error = function(e) {
          warning(paste0("error：process ", cell, "-", dataset, "-", gene, size, "fail：", e$message))
        })
      }
    }
  }
}
message("=== finish all exp similar calculate ===")





