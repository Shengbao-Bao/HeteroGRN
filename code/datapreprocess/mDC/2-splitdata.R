library(dplyr)
set.seed(123)
process_data_structured <- function(tf_file, target_file, net_file, output_prefix, dataset_type = "Specific") {
  base_dir <- "/home/bsb/mywork/data/mDC/"
  dir_map <- list(
    Specific = "SpecificDataset/groundtruth/",
    NonSpecific = "NonSpecificDataset/groundtruth/",
    String = "StringDataset/groundtruth/"  
  )
  if (!dataset_type %in% names(dir_map)) {
    stop("dataset_type must be one of: 'Specific', 'NonSpecific', 'String'")
  }
  setwd(file.path(base_dir, dir_map[[dataset_type]]))
  tf <- read.table(tf_file)[,1]
  target <- read.table(target_file)[,1]
  net <- read.table(net_file)[, 1:2, drop = FALSE]
  colnames(net) <- c("TF", "Target")
  net$Label <- 1  
  all_pairs <- expand.grid(TF = tf, Target = target)
  neg_samples <- anti_join(all_pairs, net, by = c("TF", "Target")) %>%
    mutate(
      TF = as.character(TF), 
      Target = as.character(Target)
    ) %>%
    filter(TF != Target)  
  neg_samples$Label <- 0 
  
  n_fold <- 5
  pos_folds <- split(net, sample(rep(1:n_fold, length.out = nrow(net))))
  neg_folds <- split(neg_samples, sample(rep(1:n_fold, length.out = nrow(neg_samples))))
  
  for (i in 1:n_fold) {
    val_fold <- ifelse(i == n_fold, 1, i + 1)
    test_set <- bind_rows(
      pos_folds[[i]],        
      neg_folds[[i]]         
    ) %>% sample_frac(1)     
    
    val_set <- bind_rows(
      pos_folds[[val_fold]],  
      neg_folds[[val_fold]]   
    ) %>% sample_frac(1)      
    
    train_folds <- setdiff(1:n_fold, c(i, val_fold))
    train_pos_all <- bind_rows(pos_folds[train_folds])
    train_neg_all <- bind_rows(neg_folds[train_folds])
    #
    target_size <- min(nrow(train_pos_all), nrow(train_neg_all)) 
    train_pos <- sample_n(train_pos_all, target_size)  
    train_neg <- sample_n(train_neg_all, target_size)  
    train_set <- bind_rows(train_pos, train_neg) %>% sample_frac(1)
    
    cat(paste("Fold", i, ":\n"))
    cat(paste("  test - posi:", sum(test_set$Label == 1), "neg:", sum(test_set$Label == 0), "\n"))
    cat(paste("  val - posi:", sum(val_set$Label == 1), "neg:", sum(val_set$Label == 0), "\n"))
    cat(paste("  train - posi:", nrow(train_pos), "neg:", nrow(train_neg), "\n\n"))
    
    write.table(train_set, paste0("train", output_prefix, "_fold_", i, ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    write.table(val_set, paste0("val", output_prefix, "_fold_", i, ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    write.table(test_set, paste0("test", output_prefix, "_fold_", i, ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    
    # 
    train_setpos <- subset(train_set, Label == 1)
    adj_matrix <- matrix(0, nrow = length(tf), ncol = length(target),
                         dimnames = list(tf, target))
    for (j in seq_len(nrow(train_setpos))) {
      tf_j <- train_setpos$TF[j]
      target_j <- train_setpos$Target[j]
      adj_matrix[tf_j, target_j] <- 1
    }
    colnames(adj_matrix)[1] <- paste("\t", colnames(adj_matrix)[1], sep="")
    write.table(adj_matrix, paste0("adj_matrix", output_prefix, "_fold_", i, ".txt"),
                row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    
    train_pairs <- with(train_set, paste(TF, Target, sep = "-"))
    test_pairs <- with(test_set, paste(TF, Target, sep = "-"))
    val_pairs <- with(val_set, paste(TF, Target, sep = "-"))
    
    cat(paste("Fold", i, "Overlap Check:\n"))
    cat(paste("  Train-Test Overlap:", length(intersect(train_pairs, test_pairs)), "\n"))
    cat(paste("  Train-Validation Overlap:", length(intersect(train_pairs, val_pairs)), "\n"))
    cat(paste("  Test-Validation Overlap:", length(intersect(test_pairs, val_pairs)), "\n\n"))
  }
}


process_data_structured("TF500.txt", "Target500.txt", "net500.txt", "500", dataset_type = "Specific")
process_data_structured("TF1000.txt", "Target1000.txt", "net1000.txt", "1000", dataset_type = "Specific")
process_data_structured("TF500.txt", "Target500.txt", "net500.txt", "500", dataset_type = "NonSpecific")
process_data_structured("TF1000.txt", "Target1000.txt", "net1000.txt", "1000", dataset_type = "NonSpecific")
process_data_structured("TF500.txt", "Target500.txt", "net500.txt", "500", dataset_type = "String")  
process_data_structured("TF1000.txt", "Target1000.txt", "net1000.txt", "1000", dataset_type = "String")
