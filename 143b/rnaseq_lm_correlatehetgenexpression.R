### AUTHOR: SONIA BOSCENCO
### AUGUST 13th 2024

## analyse 143b heteroplasmy using linear model
## to account for endogenous mutation heteroplasmy
rm(list=ls())

setwd("~/Desktop/reznik/1227_atac/analysis/")
source("../prerequisites.R")

save=T

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load heteroplasmy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
het_df                 <- as.data.frame(rbind(fread("~/Desktop/reznik/1227_atac//data/het/BErep1-1227.csv") %>%
                          mutate(cell = gsub("-1", "-1", cell)) %>% mutate(lib = "BErep1"),
                          fread("~/Desktop/reznik/1227_atac//data/het/BErep2-1227.csv") %>%
                          mutate(cell = gsub("-1", "-2", cell)) %>% mutate(lib = "BErep2"),
                          fread("~/Desktop/reznik/1227_atac//data/het/WT-1227.csv") %>%
                          mutate(cell = gsub("-1", "-3", cell)) %>% mutate(lib = "WT"))
)

het_df9369             <- as.data.frame(rbind(fread("~/Desktop/reznik/1227_atac//data/het/BErep1-9369.csv") %>%
                                                mutate(cell = gsub("-1", "-1", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/BErep2-9369.csv") %>%
                                                mutate(cell = gsub("-1", "-2", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/WT-9369.csv") %>%
                                                mutate(cell = gsub("-1", "-3", cell)))
)

het_df1230             <- as.data.frame(rbind(fread("~/Desktop/reznik/1227_atac//data/het/BErep1-1230.csv") %>%
                                                mutate(cell = gsub("-1", "-1", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/BErep2-1230.csv") %>%
                                                mutate(cell = gsub("-1", "-2", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/WT-1230.csv") %>%
                                                mutate(cell = gsub("-1", "-3", cell)))
)

het_df1233             <- as.data.frame(rbind(fread("~/Desktop/reznik/1227_atac//data/het/BErep1-1233.csv") %>%
                                                mutate(cell = gsub("-1", "-1", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/BErep2-1233.csv") %>%
                                                mutate(cell = gsub("-1", "-2", cell)),
                                              fread("~/Desktop/reznik/1227_atac//data/het/WT-1233.csv") %>%
                                                mutate(cell = gsub("-1", "-3", cell)))
)

write.csv(het_df, file = "~/Desktop/reznik/rrna-hotspots/data/processed/SI_figures/heteroplasmy_1227_all.csv", row.names = FALSE)

merged_het_df          <- merge(het_df, het_df9369, by = "cell")
merged_het_df          <- merge(merged_het_df, het_df1230, by = "cell")
merged_het_df          <- merge(merged_het_df, het_df1233, by = "cell")

rownames(merged_het_df)<- merged_het_df$cell
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get stats on mutant libaries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
het_df_lib1               <- het_df[het_df$lib == "BErep1", ]
het_df_lib1               <- het_df_lib1[het_df_lib1$cov >= 25, ]
het_df_lib1$af1227        <- as.numeric(het_df_lib1$af1227)

het_df_lib2               <- het_df[het_df$lib == "BErep2", ]
het_df_lib2               <- het_df_lib2[het_df_lib2$cov >= 25, ]
het_df_lib2$af1227        <- as.numeric(het_df_lib2$af1227)

median_lib1               <- median(het_df_lib1$af1227)
median_lib2               <- median(het_df_lib2$af1227)

mean(het_df_lib1$af1227)
mean(het_df_lib2$af1227)
perc_cells_less10_lib1    <- sum(het_df_lib1$af1227 < 0.1) / nrow(het_df_lib1)
perc_cells_less10_lib2    <- sum(het_df_lib2$af1227 < 0.1) / nrow(het_df_lib2)

perc_cells_greater90_lib1 <- sum(het_df_lib1$af1227 > 0.80) / nrow(het_df_lib1)
perc_cells_greater90_lib2 <- sum(het_df_lib2$af1227 > 0.80) / nrow(het_df_lib2)

cat(paste0("################## BEREP1 ##################"))
cat(paste0("Number of cells with >= 25 reads at m.1227: ", nrow(het_df_lib1)))
cat(paste0("Median Heteroplasmy: ", median_lib1))
cat(paste0("Percentage of cells with < 10% heteroplasmy: ", perc_cells_less10_lib1))
cat(paste0("Percentage of cells with > 90% heteroplasmy: ", perc_cells_greater90_lib1))

cat(paste0("################## BEREP2 ##################"))
cat(paste0("Number of cells with >= 25 reads at m.1227: ", nrow(het_df_lib2)))
cat(paste0("Median Heteroplasmy: ", median_lib2))
cat(paste0("Percentage of cells with < 10% heteroplasmy: ", perc_cells_less10_lib2))
cat(paste0("Percentage of cells with > 90% heteroplasmy: ", perc_cells_greater90_lib2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create and process seurat object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create Seurat object
gex_mat                    <- Read10X_h5("~/Desktop/reznik/1227_atac//data/aggr/filtered_feature_bc_matrix.h5")$`Gene Expression`
so                         <- CreateSeuratObject(counts = gex_mat, meta.data = merged_het_df)
so[["percent.mt"]]         <- PercentageFeatureSet(so, pattern = "^MT-")

# check summary stats
summary(so[["percent.mt"]])
summary(so[["cov1227"]])
summary(so[["nFeature_RNA"]])
VlnPlot(so, features = c("percent.mt", "cov1227", "nFeature_RNA"))

# filter seurat object
so                      <- subset(so, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 5 &
                                  cov1227 >= 25 & cov1227 < 250)
dim(so)

# Do standard dim reduction
so                    <- NormalizeData(so) %>%
                         FindVariableFeatures() %>% 
                         ScaleData() %>% 
                         RunPCA() %>% 
                         RunUMAP(dims = 1:30) %>% 
                         FindNeighbors() %>% 
                         FindClusters(resolution = 0.1)
FeaturePlot(so, "af1227", split.by = "lib")

saveRDS(so, file = "~/Desktop/reznik/rrna-hotspots/revisions/data/seurat_143b_processed.rds")

all.genes            <- rownames(gex_mat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to correlate heteroplasmy to gene expression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gene_expression_mat      <- as.data.frame(t(so@assays$RNA$data))
het                      <- as.data.frame(so$af1227)
cov                      <- as.data.frame(so$cov1227)
cluster                  <- as.data.frame(so$seurat_clusters)
lib                      <- as.data.frame(so$lib)
het9369                  <- as.data.frame(so$af9369)
het1230                  <- as.data.frame(so$af1230)
het1233                  <- as.data.frame(so$af1233)
cov9369                  <- as.data.frame(so$cov9369)
cov1230                  <- as.data.frame(so$cov1230)
cov1233                  <- as.data.frame(so$cov1233)

colnames(het)            <- c("het")
colnames(het9369)        <- c("het9369")
colnames(het1230)        <- c("het1230")
colnames(het1233)        <- c("het1233")
colnames(cluster)        <- c("cluster")
colnames(lib)            <- c("lib")
colnames(cov)            <- c("cov")
colnames(cov9369)        <- c("cov9369")
colnames(cov1230)        <- c("cov1230")
colnames(cov1233)        <- c("cov1233")

gene_expression_mat$rows <- rownames(gene_expression_mat)
het$rows                 <- rownames(het)
cluster$rows             <- rownames(cluster)
lib$rows                 <- rownames(lib)
cov$rows                 <- rownames(cov)
het9369$rows             <- rownames(het9369)
het1230$rows             <- rownames(het1230)
het1233$rows             <- rownames(het1233)
cov9369$rows             <- rownames(cov9369)
cov1230$rows             <- rownames(cov1230)
cov1233$rows             <- rownames(cov1233)

master_df                <- merge(gene_expression_mat, het, by = "rows")
master_df                <- merge(master_df, cluster, by = "rows")
master_df                <- merge(master_df, lib, by = "rows")
master_df                <- merge(master_df, cov, by = "rows")
master_df                <- merge(master_df, het9369, by = "rows")
master_df                <- merge(master_df, het1230, by = "rows")
master_df                <- merge(master_df, het1233, by = "rows")
master_df                <- merge(master_df, cov9369, by = "rows")
master_df                <- merge(master_df, cov1230, by = "rows")
master_df                <- merge(master_df, cov1233, by = "rows")

BE1                      <- subset(master_df, lib == "BErep1")
BE2                      <- subset(master_df, lib == "BErep2")
WT                       <- subset(master_df, lib == "WT")

correlate_lm_genes       <- function(gex, filename){
  pb <- txtProgressBar(min = 0, max = length(all.genes), style = 3)
  coeffs                 <- numeric((length(all.genes)))
  p_values               <- numeric((length(all.genes)))
  
  coeffs9369             <- numeric((length(all.genes)))
  p_values9369           <- numeric((length(all.genes)))
  
  names(coeffs)          <- all.genes
  names(p_values)        <- all.genes
  
  names(coeffs9369)          <- all.genes
  names(p_values9369)        <- all.genes
  
  i <- 0
  for(gene in all.genes){
    model          <- lm(gex[[gene]]~ het + het9369 + het1230 + het1233, data = gex)
    tidy_model     <- tidy(model)
    coeffs[gene]   <- tidy_model$estimate[tidy_model$term == "het"]
    coeffs9369[gene]   <- tidy_model$estimate[tidy_model$term == "het9369"]
    p_values[gene]  <- tidy_model$p.value[tidy_model$term == "het"]
    p_values9369[gene]  <- tidy_model$p.value[tidy_model$term == "het9369"]
    
    setTxtProgressBar(pb, i)
    i <- i+1
  }
  
  results_df <- data.frame(
    Gene = names(coeffs),
    Coefficient_1227 = coeffs,
    P_Value_1227 = p_values,
    Coefficient_9369 = coeffs9369,
    P_Value_9369 = p_values9369,
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df[order(results_df$Coefficient_1227, decreasing = TRUE), ]
  write.csv(results_df, file = filename, row.names = F)
  close(pb)
  return(results_df)
}

BE1_sub             <- BE1 %>% subset(cov > 10 & cov1230 > 10 & cov1233 > 10 & cov9369 > 10) 
BE2_sub             <- BE2 %>% subset(cov > 10 & cov1230 > 10 & cov1233 > 10 & cov9369 > 10) 

BE1_lm_ranks        <- correlate_lm_genes(BE1_sub, "~/Desktop/reznik/1227_atac/results/BE1_lm_ranks_NEW.csv")
BE2_lm_ranks        <- correlate_lm_genes(BE2_sub, "~/Desktop/reznik/1227_atac/results/BE2_lm_ranks.csv")

BE1_lm_coeffs <- BE1_lm_ranks$Coefficient_1227
names(BE1_lm_coeffs) <- BE1_lm_ranks$Gene
BE1_lm_coeffs <- sort(BE1_lm_coeffs, decreasing = TRUE)

BE2_lm_coeffs <- BE2_lm_ranks$Coefficient_1227
names(BE2_lm_coeffs) <- BE2_lm_ranks$Gene
BE2_lm_coeffs <- sort(BE2_lm_coeffs, decreasing = TRUE)

BE1_lm_sorted_ranks <-  sort(BE1_lm_ranks$Coefficient_1227, decreasing = TRUE)
sub_BE1 <- BE1 %>% subset(cov > 10 & cov1230 > 10 & cov1230 > 10 & cov9369 > 10 & het1230 < 0.05 & het1233 < 0.05 & het9369 < 0.05 & het > 0.01)
sub_BE2 <- BE2 %>% subset(cov > 10 & cov1230 > 10 & cov1230 > 10 & cov9369 > 10 & het1230 < 0.05 & het1233 < 0.05 & het9369 < 0.05 & het > 0.01)


# gsea 

gmtfile                  <- "~/Desktop/reznik/reference/gsea_terms/c5.go.v2023.1.Hs.symbols.gmt"
#gmtfile                  <- "~/Desktop/reznik/reference/gsea_terms/h.all.v7.2.symbols.gmt"
gs                       <-  read.gmt(gmtfile)
gsea                   <- function(gs, ranks){
  resgsea = fgsea(pathways = gs,stats = ranks,minSize=1,maxSize=500)
  resgsea = resgsea[order(resgsea$padj,decreasing = FALSE),]
  resdf = apply(resgsea,2,as.character)
  resdf[,'padj'] <- as.numeric(resdf[,'padj'])
  resdf[,'NES'] <- as.numeric(resdf[,'NES'])
  #resdf_sig = resdf[which(resdf[,'padj'] < .05),]
  return(resdf)
}


BE1_lm_gsea <- gsea(gs, BE1_lm_coeffs) 
BE2_lm_gsea <- gsea(gs, BE2_lm_coeffs) 
BE1_lm_gsea <- as.data.frame(BE1_lm_gsea)
BE2_lm_gsea <- as.data.frame(BE2_lm_gsea)

BE1_lm_gsea$padj <- as.numeric(BE1_lm_gsea$padj)
BE2_lm_gsea$padj <- as.numeric(BE2_lm_gsea$padj)
