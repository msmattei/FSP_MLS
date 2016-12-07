
## Tissue specificity score calculated for each experiment (GSE..)
## for future analysis only the tissue specificity score calculated with the experiment with more tissues are used
## the better way to calculate the tissue specificity score using different experiment should be tested in the future

directory <- "C://Mimi/UNI/Master/MLS_BIOINFORMATICS/Isemestre/First_step_project_MLS/"
setwd(directory)


# Tissue specific score ----------------------------------------

tau <- function(rpkm){
  score <- sum(1 - (rpkm / max(rpkm))) / (length(rpkm) - 1)
  return(score)
}


# Compute -----------------------------------------------------------------

## for human

experiment_tissue <- list.files("Data/Human/") # list of all GSE..tissues.RData created with data_extraction.R

for (gse in experiment_tissue){
  load(paste0("Data/Human/", gse))
  paralog_gene_ID <- unique(tissues[[1]]$Ensembl_Gene_ID_short)
  
  tissue_spec_score <- NULL
  for (gene in paralog_gene_ID){
    gene_expression <- NULL
    for (tiss in names(tissues)){
      rpkm <- mean(tissues[[tiss]]$RPKM_short[tissues[[tiss]]$Ensembl_Gene_ID_short==gene])
      gene_expression <- c(gene_expression,rpkm)
    }
    score <- tau(gene_expression)
    tissue_spec_score <- c(tissue_spec_score, score)
  }
  
  tissue_spec_score <- data.frame(Gene.ID = paralog_gene_ID, tissue_spec_score)[complete.cases(data.frame(Gene.ID = paralog_gene_ID, tissue_spec_score)),]
  write.table(tissue_spec_score, paste0("Output/Human/", strsplit(gse, "[.]")[[1]][1], "spec_score", Sys.Date(), ".txt"), col.names = F, row.names = F, quote = F)
  
}


## for mouse

mouse_tissue <- list.files("Data/Mouse/")

for (gse in mouse_tissue){
  load(paste0("Data/Mouse/", gse))
  paralog_gene_ID <- unique(tissues[[1]]$Ensembl_Gene_ID_short)
  
  tissue_spec_score <- NULL
  for (gene in paralog_gene_ID){
    gene_expression <- NULL
    for (tiss in names(tissues)){
      rpkm <- mean(tissues[[tiss]]$RPKM_short[tissues[[tiss]]$Ensembl_Gene_ID_short==gene])
      gene_expression <- c(gene_expression,rpkm)
    }
    score <- tau(gene_expression)
    tissue_spec_score <- c(tissue_spec_score, score)
  }
  
  tissue_spec_score <- data.frame(Gene.ID = paralog_gene_ID, tissue_spec_score)[complete.cases(data.frame(Gene.ID = paralog_gene_ID, tissue_spec_score)),]
  write.table(tissue_spec_score, paste0("Output/Mouse/", strsplit(gse, "[.]")[[1]][1], "spec_score", Sys.Date(), ".txt"), col.names = F, row.names = F, quote = F)
  
}
