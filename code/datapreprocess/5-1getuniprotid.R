library(biomaRt)
#read human tf
hESC_spec  <- read.table("/home/bsb/mywork/data/hESC/SpecificDataset/groundtruth/TF1000.txt")[,1]
hESC_nonspec <- read.table("/home/bsb/mywork/data/hESC/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
hESC_string  <- read.table("/home/bsb/mywork/data/hESC/StringDataset/groundtruth/TF1000.txt")[,1]
hHep_spec    <- read.table("/home/bsb/mywork/data/hHep/SpecificDataset/groundtruth/TF1000.txt")[,1]
hHep_nonspec <- read.table("/home/bsb/mywork/data/hHep/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
hHep_string  <- read.table("/home/bsb/mywork/data/hHep/StringDataset/groundtruth/TF1000.txt")[,1]
human_tf_all <- unique(c(hESC_spec, hESC_nonspec, hESC_string,
                         hHep_spec, hHep_nonspec, hHep_string))
# 1.2 BioMart  (Swiss-Prot)
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_bm <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  filters = "hgnc_symbol",
  values = human_tf_all,
  mart = ensembl_human
)
human_bm <- human_bm[human_bm$uniprotswissprot != "", ]
# 1.3 supply （ TSV.GZ）
human_tsv <- read.delim("/home/bsb/mywork/rawdata/structure/idmapping_model_organism_9606_AND_revie_2026_01_06.tsv.gz",
                        header = TRUE, stringsAsFactors = FALSE)
colnames(human_tsv)[1:2] <- c("hgnc_symbol", "uniprotswissprot")
# 1.4 merge
human_uniprot_all <- rbind(human_bm, human_tsv[, 1:2])
# 1.5 Filter: Keep only UniProt IDs with structural files
human_structure_dir <- "/home/bsb/mywork/rawdata/structure/human/"
human_files <- list.files(human_structure_dir)
human_uniprot_in_files <- sapply(strsplit(human_files, "-"), 
                                 function(x) if (length(x) >= 2) x[2] else NA)
human_uniprot_final <- human_uniprot_all[human_uniprot_all$uniprotswissprot %in% human_uniprot_in_files, ]
setdiff(human_tf_all, human_uniprot_final$hgnc_symbol)
# 1.6 save result
write.table(human_uniprot_final, 
            "/home/bsb/mywork/data/middlefile/human_tf_uniprot.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(human_uniprot_final$uniprotswissprot,
            "/home/bsb/mywork/data/middlefile/human_uniprot_ids.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



# ===================================================================
# ===================================================================
# 2. (mouse)
# 2.1 read mouse tf list
mESC_spec  <- read.table("/home/bsb/mywork/data/mESC/SpecificDataset/groundtruth/TF1000.txt")[,1]
mESC_nonspec <- read.table("/home/bsb/mywork/data/mESC/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
mESC_string  <- read.table("/home/bsb/mywork/data/mESC/StringDataset/groundtruth/TF1000.txt")[,1]
mDC_spec     <- read.table("/home/bsb/mywork/data/mDC/SpecificDataset/groundtruth/TF1000.txt")[,1]
mDC_nonspec  <- read.table("/home/bsb/mywork/data/mDC/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
mDC_string   <- read.table("/home/bsb/mywork/data/mDC/StringDataset/groundtruth/TF1000.txt")[,1]
mHSC_L_spec  <- read.table("/home/bsb/mywork/data/mHSC-L/SpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_L_nonspec <- read.table("/home/bsb/mywork/data/mHSC-L/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_L_string  <- read.table("/home/bsb/mywork/data/mHSC-L/StringDataset/groundtruth/TF1000.txt")[,1]
mHSC_GM_spec   <- read.table("/home/bsb/mywork/data/mHSC-GM/SpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_GM_nonspec <- read.table("/home/bsb/mywork/data/mHSC-GM/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_GM_string  <- read.table("/home/bsb/mywork/data/mHSC-GM/StringDataset/groundtruth/TF1000.txt")[,1]
mHSC_E_spec    <- read.table("/home/bsb/mywork/data/mHSC-E/SpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_E_nonspec <- read.table("/home/bsb/mywork/data/mHSC-E/NonSpecificDataset/groundtruth/TF1000.txt")[,1]
mHSC_E_string  <- read.table("/home/bsb/mywork/data/mHSC-E/StringDataset/groundtruth/TF1000.txt")[,1]
mouse_tf_all <- unique(c(mESC_spec, mESC_nonspec, mESC_string,
                         mDC_spec, mDC_nonspec, mDC_string,
                         mHSC_L_spec, mHSC_L_nonspec, mHSC_L_string,
                         mHSC_GM_spec, mHSC_GM_nonspec, mHSC_GM_string,
                         mHSC_E_spec, mHSC_E_nonspec, mHSC_E_string))

# 2.2 BioMart (Swiss-Prot)
ensembl_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_bm <- getBM(
  attributes = c("mgi_symbol", "uniprotswissprot"),
  filters = "mgi_symbol",
  values = mouse_tf_all,
  mart = ensembl_mouse
)
mouse_bm$mgi_symbol <- toupper(mouse_bm$mgi_symbol)
mouse_bm <- mouse_bm[mouse_bm$uniprotswissprot != "", ]
# 2.3 supply mapping（TSV.GZ）
mouse_tsv <- read.delim("/home/bsb/mywork/rawdata/structure/idmapping_Mus_musculus_Mouse_2026_01_06.tsv.gz",
                        header = TRUE, stringsAsFactors = FALSE)
mouse_tsv <- subset(mouse_tsv, Organism == "Mus musculus (Mouse)")
colnames(mouse_tsv)[1:2] <- c("mgi_symbol", "uniprotswissprot")
# 2.4 merge all mapping
mouse_uniprot_all <- rbind(mouse_bm, mouse_tsv[, 1:2])
# 2.5 first epoch：preserve with structure 
mouse_structure_dir <- "/home/bsb/mywork/rawdata/structure/mouse/"
mouse_files <- list.files(mouse_structure_dir)
mouse_uniprot_in_files <- sapply(strsplit(mouse_files, "-"), 
                                 function(x) if (length(x) >= 2) x[2] else NA)
mouse_with_structure <- mouse_uniprot_all[mouse_uniprot_all$uniprotswissprot %in% mouse_uniprot_in_files, ]
# 2.6 second：supply missing tf
missing_tf <- setdiff(mouse_tf_all, mouse_with_structure$mgi_symbol)
# 2.61 : BioMart  Swiss-Prot
bm_missing <- mouse_bm[mouse_bm$mgi_symbol %in% missing_tf, ]
# 2.62: TSV  Reviewed 
tsv_reviewed <- subset(mouse_tsv,mgi_symbol %in% missing_tf & Reviewed == "reviewed")
# 2.63
remaining_missing <- setdiff(missing_tf, c(bm_missing$mgi_symbol, tsv_reviewed$mgi_symbol))
tsv_unreviewed <- mouse_tsv[mouse_tsv$mgi_symbol %in% remaining_missing, ]
tsv_unreviewed_first <- tsv_unreviewed[!duplicated(tsv_unreviewed$mgi_symbol), ]
# 2.7 merge , final result
mouse_uniprot_final <- rbind(
  mouse_with_structure,
  bm_missing[, 1:2],
  tsv_reviewed[, 1:2],
  tsv_unreviewed_first[, 1:2]
)
setdiff(mouse_tf_all, mouse_uniprot_final$mgi_symbol)
sort(table(mouse_uniprot_final$mgi_symbol)) 
#HBP1:Q8R316 Q9D1N2(Wrong)
#MLL2:O08550 Q6PDK2(mll4)
#PMS1:Q5NC83 Q8K119(unreview)
mouse_uniprot_final <- mouse_uniprot_final[
  !mouse_uniprot_final$uniprotswissprot %in% c("Q9D1N2", "Q6PDK2", "Q8K119"),
]
# 2.8 save
write.table(mouse_uniprot_final, 
            "/home/bsb/mywork/data/middlefile/mouse_tf_uniprot.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(mouse_uniprot_final$uniprotswissprot,
            "/home/bsb/mywork/data/middlefile/mouse_uniprot_ids.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
missing_structre=setdiff(mouse_uniprot_final$uniprotswissprot,mouse_with_structure$uniprotswissprot)
write.table(missing_structre,
            "/home/bsb/mywork/data/middlefile/missing_struture_ids.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)














