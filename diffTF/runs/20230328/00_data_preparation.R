# preparing the sampleData.tsv and an RNA-seq counts table.

#loading the RData file from ab2_nk_tumor_project 
#maybe some libraries will be missing

install.packages("DESeq2")
install.packages("stringr")
library(DESeq2)
load("/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab2_20230206_nk_tumor_project/old_code/NK_Tum_immunoedit/results_dir/nk_tum_immunoedit_complete_ab2_AB_CD_sep//nk_tum_immunoedit_complete_ab2_AB_CD_sep_dds_objects.RData")

#get the counts table for diffTFs  AB cell lines 
fl_smpl_nm <- gsub("_Day\\d\\d","",dds_subset@colData@listData[["full_sample_name"]])
fl_smpl_nm <- gsub("_\\d_","_",fl_smpl_nm)

sample_name <- paste0(dds_subset@colData@listData[["original_experiment"]], "_",
                      fl_smpl_nm, "_",
                      dds_subset@colData@listData[["day_cell_harvesting"]], "_REP",
                      dds_subset@colData@listData[["technical_replicate"]])

sample_name <- gsub("^EXP\\s9.", "EXP9_", sample_name)
sample_name <- gsub("\\+", "_plus_", sample_name)
sample_name <- gsub("_Tumor_only_Day\\d_", "_Tumor_only_", sample_name)

#These samples from RNA-seq are missing in ATAC data.
#EXP9_6_13G_Tumor_plus_NK_Day14_REP1 is missing in ./
#EXP9_6_13G_Tumor_plus_NK_Day14_REP2 is missing in ./
#EXP9_6_13G_Tumor_plus_NK_Day14_REP3 is missing in ./
#EXP9_6_13H_Tumor_plus_NK_Day14_REP1 is missing in ./
#EXP9_6_13H_Tumor_plus_NK_Day14_REP2 is missing in ./
#EXP9_6_13H_Tumor_plus_NK_Day14_REP3 is missing in ./
#
#Instead there is Day_10 version.
#
#I'm going to replace in the RNA-seq count table the day 14 with day 10, so that samples match.

sample_name <- gsub("EXP9_6_13G_Tumor_plus_NK_Day14", "EXP9_6_13G_Tumor_plus_NK_Day10", sample_name)
sample_name <- gsub("EXP9_6_13H_Tumor_plus_NK_Day14", "EXP9_6_13H_Tumor_plus_NK_Day10", sample_name)
col_IDS <- c("ENSEMBL", sample_name)
RNA_seq_raw_counts <- dds_subset@assays@data@listData[["counts"]]
RNA_seq_raw_counts <- cbind(row.names(RNA_seq_raw_counts), RNA_seq_raw_counts)
RNA_seq_raw_counts <- rbind(col_IDS,RNA_seq_raw_counts)

write.table(RNA_seq_raw_counts, file = "~/workspace/RNA_seq_raw_counts_AB.tsv",
            sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Building the sample table
sampleID_list <- list("sampleID" = sample_name)
bamReads_list <- list("bamReads" = paste0("/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab3_20230220_nk_tumor_project/data_ATACSeq_nfcore_run/results/bwa/merged_library/", sample_name, ".mLb.clN.sorted.bam"))
cond_summ <- sample_name
cond_summ[stringr::str_detect(cond_summ, pattern = "_Tumor_only_")] <- "Tumor_only_TP_1"
cond_summ[stringr::str_detect(cond_summ, pattern = "_Tumor_plus_NK")] <- "Tumor_pls_NK_TP2"
conditionSummary_list <- list("conditionSummary" = cond_summ)
experiment_id_list <- list("experiment_id" = str_extract(sample_name, "EXP\\d_\\d"))

sampleData.tsv <- data.frame(sampleID_list, 
                             bamReads_list, 
                             conditionSummary_list,
                             experiment_id_list)
write.table(sampleData.tsv, file = "./sampleData.tsv",
            sep ="\t", quote = FALSE, row.names = FALSE)





#### Preparing consensus peaks table ####
dds_consens_1 <- readRDS(file = "~/workspace/Rproject/Results/RDSs/ATACSeq_DESeq2_obj_complete_consensus_1.rds")
dds_cons_1_df <- data.frame(dds_consens_1@rowRanges@elementMetadata@listData[["Chr"]],
                            dds_consens_1@rowRanges@elementMetadata@listData[["Start"]],
                            dds_consens_1@rowRanges@elementMetadata@listData[["End"]],
                            dds_consens_1@rowRanges@elementMetadata@listData[["PeakId"]])
write.table(dds_cons_1_df, file = "~/workspace/references/dds_consensus_1_peaks.tsv",
            sep = '\t', col.names = FALSE, quote = FALSE, row.names = FALSE)
