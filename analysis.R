# Load Packages
library(biomaRt)
library(DESeq2)
library(tidyverse)
library(ashr)
library(data.table)


## Determine DEGs for duck and chicken

# Filter count and experimental data for tissues
returnPart <- function(data, exp, part){
  exp_part <- exp %>% filter(Factor.Value.organism.part. == part) %>%
    mutate(combinedFactor = as.factor(gsub(" ", "",
                                           paste(Factor.Value.infect.,
                                                  Factor.Value.time., sep = '_'))))
  data_part <- data %>% select(rownames(exp_part))
  exp_part = exp_part %>% slice(match(colnames(data_part), rownames(exp_part)))
  return(list("data" = data_part, "exp" = exp_part))
}

# Run DESeq for both tissues
run_DESeq <- function(counts, exp){
  lung <- returnPart(counts, exp, "lung")
  ileum <- returnPart(counts, exp, "ileum")
  dds_lung <- DESeqDataSetFromMatrix(countData = lung$data,
                                          colData = lung$exp, design = ~ combinedFactor) %>%
    DESeq()
  dds_ileum <- DESeqDataSetFromMatrix(countData = ileum$data,
                                           colData = ileum$exp, design = ~ combinedFactor) %>%
    DESeq()
  return(list(lung = dds_lung, ileum = dds_ileum))
}

# Retrieve results from different contrasts
get_results <- function(dds, contrast, type = "normal", padj.cutoff = 0.05,
                        lfc.cutoff = 0.58, alpha = 0.05){
  res <- results(dds, contrast = contrast, alpha = alpha)
  res = lfcShrink(dds, contrast = contrast, res = res, type = type)
  res_df <- data.frame(res)
  res_df <- rownames_to_column(res_df, var = "gene")

  return(res_df)
}

# Load data
duck_exp <- read.table("Duck/Experiment/E-MTAB-2909-experiment-design.tsv",
                       sep = "\t", header = TRUE, row.names = 1)
duck_counts <- read.table("Duck/Experiment/E-MTAB-2909-raw-counts.tsv",
                          sep = "\t", header = TRUE, row.names=1)  %>% select(-1)
chicken_exp <- read.table("Chicken/Experiment/E-MTAB-2908-experiment-design.tsv",
                          sep = "\t", header = TRUE, row.names = 1)
chicken_counts <- read.table("Chicken/Experiment/E-MTAB-2908-raw-counts.tsv",
                             sep = "\t", header = TRUE, row.names = 1) %>% select(-1)

# Define contrasts
contrast <- list()
contrast[["H5N1_ileum_1day"]] <- c("combinedFactor", "H5N1_1day", "none_1day")
contrast[["H5N2_ileum_1day"]] <- c("combinedFactor", "H5N2_1day", "none_1day")
contrast[["H5N1_ileum_3day"]] <- c("combinedFactor", "H5N1_3day", "none_3day")
contrast[["H5N2_ileum_3day"]] <- c("combinedFactor", "H5N2_3day", "none_3day")
contrast[["H5N1_lung_1day"]] <- c("combinedFactor", "H5N1_1day", "none_1day")
contrast[["H5N2_lung_1day"]] <- c("combinedFactor", "H5N2_1day", "none_1day")
contrast[["H5N1_lung_3day"]] <- c("combinedFactor", "H5N1_3day", "none_3day")
contrast[["H5N2_lung_3day"]] <- c("combinedFactor", "H5N2_3day", "none_3day")

padj.cutoff <- 0.05
lfc.cutoff <- 0.58
type = "ashr"

# Run DESeq for chicken and duck
duck_dds <- run_DESeq(duck_counts, duck_exp)
chicken_dds <- run_DESeq(chicken_counts, chicken_exp)

# Retrieve results for all contrasts and write to storage
duck_res <- list()
duck_degs <- list()
chicken_res <- list()
chicken_degs <- list()

for(contr in names(contrast)){
  #ileum else lung
  if(grepl('ileum', contr)){
    duck_res[[contr]] <- get_results(dds = duck_dds$ileum, contrast = contrast[[contr]], type = type)
    duck_degs[[contr]] <- filter(duck_res[[contr]], padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
    chicken_res[[contr]] <- get_results(dds = chicken_dds$ileum, contrast = contrast[[contr]], type = type)
    chicken_degs[[contr]] <- filter(chicken_res[[contr]], padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
  }else{
    duck_res[[contr]] <- get_results(dds = duck_dds$lung, contrast = contrast[[contr]], type = type)
    duck_degs[[contr]] <- filter(duck_res[[contr]], padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
    chicken_res[[contr]] <- get_results(dds = chicken_dds$lung, contrast = contrast[[contr]], type = type)
    chicken_degs[[contr]] <- filter(chicken_res[[contr]], padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
  }
  file_name <- paste0('Chicken/DEG/', contr, '.txt')
  fwrite(chicken_degs[[contr]], file_name, sep = '\t', col.names = TRUE)
  file_name <- paste0('Duck/DEG/', contr, '.txt')
  fwrite(duck_degs[[contr]], file_name, sep = '\t', col.names = TRUE)
}



## Determine orthologous genes by querying BioMart
# To reproduce the results of the paper, choose the corresponding 
# ensembl release from the archive (CAU_duck1.0, GRCg6a)
duck_mart <- useMart("ensembl", dataset="applatyrhynchos_gene_ensembl")
chicken_mart <- useMart("ensembl", dataset="ggallus_gene_ensembl")
attr <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")

# Retrieve orthologous genes for every contrast and write to storage
chicken_degs_in_duck <- list()
duck_degs_in_chicken <- list()

for(contr in names(contrast)){
  if(nrow(chicken_degs[[contr]]) != 0){
    chicken_degs_in_duck[[contr]] <- getLDS(attributes = "ensembl_gene_id",
                                            filters = "ensembl_gene_id",
                                            values = chicken_degs[[contr]]$gene, mart = chicken_mart,
                                            attributesL = "ensembl_gene_id", martL = duck_mart)
    colnames(chicken_degs_in_duck[[contr]])  <- c("chicken_genes", "duck_genes")
    write.table(subset(chicken_degs_in_duck[[contr]], !duplicated(duck_genes)),
                paste0("Chicken/DEG_orth/", contr),
                col.names = T, sep = "\t", row.names = F, quote = F)
  }else{
    write.table(data.frame(), paste0("Chicken/DEG_orth/", contr),
                col.names = T, sep = "\t", row.names = F, quote = F)
  }
  
  duck_degs_in_chicken[[contr]] <- getLDS(attributes = "ensembl_gene_id", mart = chicken_mart,
                                 attributesL = "ensembl_gene_id", filtersL = "ensembl_gene_id",
                                 valuesL = duck_degs[[contr]]$gene, martL = duck_mart)
  
  colnames(duck_degs_in_chicken[[contr]])  <- c("chicken_genes", "duck_genes")
  write.table(subset(duck_degs_in_chicken[[contr]], !duplicated(chicken_genes)),
              paste0("Duck/DEG_orth/", contr),
              col.names = T, sep = "\t", row.names = F, quote = F)
}


## Prepare data for TF enrichment
# Retrieve all genes for chicken and duck from BioMart
chicken_genes <- getBM(attributes = attr, mart = chicken_mart)
duck_genes <- getBM(attributes = attr, mart = duck_mart)
# Consider only H5N1
relevant_degs_duck <- duck_degs_in_chicken[c(1,3,5,7)]

# Create gene list for duck DEGs and orthologous genes in chicken
for(name in names(relevant_degs_duck)){
  duck_degs_list <- filter(duck_genes, ensembl_gene_id %in% relevant_degs_duck[[name]]$duck_genes)
  chicken_orthologues_degs_list <- filter(chicken_genes, ensembl_gene_id %in% relevant_degs_duck[[name]]$chicken_genes)
  
  write.table(duck_degs_list, paste0("Duck/Promoter/Info/", name, ".txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(chicken_orthologues_degs_list, paste0("Chicken/Promoter/Info/", name, ".txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
} 

                                                                                                                                                       
# The following steps are described in README.md
# All result files are provided in the repository


