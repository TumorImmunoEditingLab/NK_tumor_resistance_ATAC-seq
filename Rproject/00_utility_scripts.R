# Utility Functions


generatePCA <- function(transf_object=NULL, cond_interest_varPart=NULL, color_variable=NULL, shape_variable=NULL, ntop_genes=500){
  
  library(DESeq2)
  library(ggplot2)
  
  #transf_object_counts <- assay(transf_object)
  #ntop_genes=nrow(transf_object_counts) 
  pcaData <- DESeq2::plotPCA(transf_object, intgroup=cond_interest_varPart, returnData=TRUE, ntop=ntop_genes) 
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=!!sym(color_variable), shape=!!sym(shape_variable))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #+ coord_fixed()
    theme_bw() 
  
}

meanExprsPerGroup <- function(dds_object=NULL, 
                              #condition_test=NULL,
                              variable=NULL,
                              cond_numerator=NULL, 
                              cond_denominator=NULL){
  require(DESeq2)
  require(dplyr)
  # extracts expression matrix (all samples) for numerator and denominator
  # calculates mean expression per group (numerator, denominator)
  # [ ] add check that there is filenames and condition column!
  # Normalized counts and means per group
  # creating only subset for actual comparison
  # may need to change condition to some other variable
  # res_extract <- condition_test
  # cond_numerator <- gsub(pattern = paste0("(", variable,"_)(.+)(_vs_.+)"), replacement = "\\2", res_extract)
  # cond_denominator <- gsub(pattern = paste0("(", variable,"_)(.+_vs_)(.+)"), replacement = "\\3", res_extract)
  
  new_sample_names <- as.data.frame(colData(dds_object)) 
  
  if("filenames" %in% names(new_sample_names)){
    new_sample_names <- new_sample_names %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
    
  } else {
    new_sample_names <- new_sample_names %>%
      tibble::rownames_to_column(., var = "filenames") %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
  }
  
  # check if variable is factor otherwise arranging in the next step will be "random"
  #if (is.factor()) {
  #  
  #}
  
  # extract filenames for each of the conditions and pivot_table
  new_sample_names <- new_sample_names %>%
    dplyr::select(filenames, !!as.name(variable)) %>% # selecting filenames and variable of interest (e.g.)
    dplyr::filter(!!as.name(variable) == cond_numerator | !!as.name(variable) == cond_denominator) %>% # filtering to keep only numerator and denominator samples
    dplyr::transmute(filenames,
                     denom_num_extract = factor(!!as.name(variable), levels=c(cond_denominator, cond_numerator))) %>% # re-factoring denominator, then numerator (but this should have been done in cond_data already)
    dplyr::arrange(., denom_num_extract) %>% # convert the strings to names with as.name and !! unquote (bang-bang) !!as.name(variable)
    # reorder variable according to numerator, denominator; so the expression output is in correct order - this assumes previous correct ordering
    dplyr::mutate(new_name = paste(denom_num_extract, filenames, sep="-")) # creating new column with variable of interest and filename
  
  
  normalized_counts <- NULL # just to make sure it does not exist from previous run
  normalized_counts <- data.frame(counts(dds_object, normalized = TRUE))
  # extracting subset
  normalized_counts <- normalized_counts %>%
    dplyr::select(new_sample_names$filenames) # keeping only conditions that are being compared and following previous order
  
  colnames(normalized_counts) <- new_sample_names$new_name
  normalized_counts <- normalized_counts %>%
    tibble::rownames_to_column("PeakId")
  
  normalized_counts_AddedMean <- normalized_counts %>%
    dplyr::mutate(., 
                  MeanExpr_denominator = rowMeans(dplyr::select(., matches(paste0(cond_denominator,"-"))), na.rm = TRUE),
                  MeanExpr_numerator = rowMeans(dplyr::select(., matches(paste0(cond_numerator,"-"))), na.rm = TRUE)) # more robust regex?!
  
  # rename mean_numerator, mean_denominator in the final column
  # use rename_all! rename(new_sample_names$new_sample_names)
  # rename colnames to condition_Sample number name
  # calculate mean across normalized counts
  
  return(normalized_counts_AddedMean)
}



extract_results_DDS <- function(dds_object = NULL,
                                coeff_name = NULL,
                                cond_numerator = NULL,
                                cond_denominator = NULL,
                                cond_variable = NULL,
                                padj_cutoff = 0.5,
                                log2FC_cutoff = 0.58){
  
  #This is the start of the PR function that I re-purposed.
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          name=coeff_name,   
                                          parallel = TRUE, 
                                          alpha = padj_cutoff)
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,   
                                 res=res_table_unshrunken, 
                                 type = "apeglm")
  row.names(res_table) <-   dds_object@rowRanges@elementMetadata@listData[["PeakId"]] #DATASET   SPECIFIC!
  
  normalized_counts_AddedMean <-   meanExprsPerGroup(dds_object=dds_object,
                                                     cond_numerator=cond_numerator,
                                                     cond_denominator=cond_denominator,
                                                     variable = cond_variable)
  normalized_counts_AddedMean$PeakId <-   dds_object@rowRanges@elementMetadata@listData[["PeakId"]]
  
  ensemblAnnot <- as.data.frame(dds_object@rowRanges@elementMetadata@listData[1:19])
  ######
  #Operations with table, merging with metadata.
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("PeakId") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, 
                                    -2^(abs(log2FoldChange)), 
                                    2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "PeakId") %>% # adding annotationl   entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "PeakId") %>% # adding   normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(PeakId, Chr, Start, End, Strand,
                  Gene.Name, Gene.Type, Annotation, Detailed.Annotation,
                  Distance.to.TSS, Nearest.PromoterID, Entrez.ID,   Nearest.Unigene,
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, pvalue, padj,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns   not mentioned; arts_with(paste0(cond_denominator,"_S") may not work if other   references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs =   funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find   better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs =   funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design =   paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes =   nrow(results_data_annot_signif),
                                        signif_genes_UP =   sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN =   sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) &   (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  dds_consens_results_list<-list(results_signif=results_data_annot_signif  , de_details=temp_results_summary_df, results_all=results_data_annot)
  return(dds_consens_results_list)
  #END OF PR FUNCTION
  
}


plotVolcano <- function(dds_results_obj=NULL, genes_of_interest=NULL,genes_of_interest2=NULL, genes_of_interest3=NULL, plot_title=NULL, max.overlaps = 10){
  results_data_annot_forPlot <- dds_results_obj %>% dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "ATAC_UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "ATAC_DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$Gene.Name %in% genes_of_interest2] <- "ATAC_UP"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$Gene.Name %in% genes_of_interest3] <- "ATAC_DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$Gene.Name %in% genes_of_interest] <- "RNAseq"
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE, levels = c("NO", "RNAseq", "ATAC_UP", "ATAC_DOWN"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point() +
    ggrepel::geom_label_repel(data = . %>% filter(Gene.Name %in% c(genes_of_interest, genes_of_interest2 ,genes_of_interest3)) %>% filter(signif_DE %in% c("RNAseq", "ATAC_UP", "ATAC_DOWN")), aes(label = Gene.Name),
                              show.legend = FALSE, max.overlaps = max.overlaps, box.padding = 0.5, segment.color = "black", min.segment.length = 0) +
    #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
    #ggrepel::geom_label_repel(data = . %>% filter(mgi_symbol %in% genes_of_interest), aes(label = mgi_symbol),
    #show.legend = FALSE) +
    #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
    geom_hline(yintercept=-log10(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
    #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
    scale_color_manual(values=c(RNAseq="orange", NO = "grey", ATAC_UP = "#DD3344", ATAC_DOWN = "#553388")) +
    theme_bw(base_size = 14) +
    # labs(x = "log2FC") + 
    ggtitle(plot_title)
  #y = "-log10( p-value )",color = "signif. DE") 
}


generatePCA_plus_shape <- function (transf_object = NULL, cond_interest_varPart = NULL, 
                                    color_variable = NULL, shape_variable = NULL,   ntop_genes = 500) 
{
  pcaData <- DESeq2::plotPCA(transf_object, intgroup = cond_interest_varPart, 
                             returnData = TRUE, ntop = ntop_genes)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData <- cbind(pcaData, list('Cell_line_label_and_exp_id' = paste(pcaData$cell_line_id, pcaData$experiment)))
  ggplot(pcaData, aes(PC1, PC2, color = !!sym(color_variable), 
                      shape = !!sym(shape_variable))) + geom_point(size = 4) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    theme_bw()
}


plotVolcano_repel <- function(dds_results_obj=NULL, genes_of_interest=NULL, plot_title=NULL){
  results_data_annot_forPlot <- dds_results_obj %>%
    dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "DOWN"
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE,
                                                 levels = c("NO", "DOWN", "UP"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point() +
    #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
    ggrepel::geom_label_repel(data = . %>% filter(mgi_symbol %in% genes_of_interest) %>% filter(signif_DE %in% c("DOWN", "UP")), aes(label = mgi_symbol), show.legend = FALSE, box.padding = 0.5, segment.color = "black", min.segment.length = 0) +
    #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
    #geom_hline(yintercept=-log10(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
    #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
    scale_color_manual(values=c(DOWN="#553388", UP="#DD3344", NO = "grey")) +
    theme_bw(base_size = 14) +
    labs(x = "log2FC") + 
    ggtitle(plot_title)
  #y = "-log10( p-value )",color = "signif. DE") 
}




#' rm_blacklisted_regions
#'
#' Removes blacklisted regions from peaks
#'
#' @param peak_gr peaks GRanges object
#' @param blacklisted_regions blacklisted regions GRanges object 
#'
#' @return GRanges object of peaks with blacklisted regions removed
#'
#' @examples rm_blacklisted_regions(peak_gr = peak_grlist[[sample_name]]
remove_blacklisted <- function(peak_gr=NULL, blacklisted_regions=blacklist){
  import::here(.from = IRanges, overlapsAny)  # IRanges::overlapsAny
  # removing regions overlapping with blacklisted
  # see https://biodatascience.github.io/compbio/bioc/ranges.html
  # If we just wanted to subset to the genes which overlap a given range, we can use overlapsAny:
  #g[overlapsAny(g, query[1])]
  #This is equivalent to the following:
  #g[g %over% query[1]] 
  
  blacklisted = sum(overlapsAny(peak_gr, blacklisted_regions))
  blacklisted_perc = (blacklisted/length(peak_gr))*100
  #not_blacklisted = sum(!overlapsAny(peak_gr, blacklisted_regions))
  message("Blacklisted regions regions: ", blacklisted, " (", round(blacklisted_perc,1), "%)")
  
  return(peak_gr[!overlapsAny(peak_gr, blacklisted_regions)])
}


#' add_flag
#'
#' Add flags to pheatmap
#' @description ADD more detailed description
#'
#' @param pheatmap refers to the original heatmap produced from the pheatmap() function
#' @param kept.labels should be a vector of labels you wish to show
#' @param repel.degree is a number in the range [0, 1], controlling how much the labels are spread out from one another
#' @param hiden.labels.character character to hide labels; e.g. if need to add spacing?!
#'
#' @return updated pheatmap
#'
#' @references original function from Z.Lin https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
#'
#' @examples add_flag(signifLog2Promoters_heatmap, kept.labels = highlight_genes, repel.degree = 0.2,hiden.labels.character = "")
add_flag <- function(pheatmap,
                     kept.labels,
                     repel.degree,
                     hiden.labels.character="") {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  # PR: hiden.labels.character - character to hide labels; e.g. if need to add spacing?!
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  # new.label$label <- ifelse(new.label$label %in% kept.labels, 
  #                           new.label$label, "")
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, hiden.labels.character)
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != hiden.labels.character)
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != hiden.labels.character],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != hiden.labels.character] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


