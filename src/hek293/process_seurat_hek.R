## AUTHOR: SONIA BOSCENCO

rm(list=ls())
library(Seurat)
setwd("~/Desktop/reznik/1227_atac/analysis/")
source("~/Desktop/reznik/1227_atac/prerequisites.R")

save=T

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load heteroplasmy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
het_df                 <- as.data.frame(rbind(fread("~/Desktop/reznik/rrna-hotspots/hek_cells/data/Experiment3_HEK293_1227G_A.csv"))
)
rownames(het_df)<- het_df$cell
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create and process seurat object
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create Seurat object
gex_mat                    <- Read10X_h5("~/Desktop/reznik/rrna-hotspots/hek_cells/data/filtered_feature_bc_matrix.h5")$`Gene Expression`
so                         <- CreateSeuratObject(counts = gex_mat, meta.data = het_df)
so[["percent.mt"]]         <- PercentageFeatureSet(so, pattern = "^MT-")

# check summary stats
summary(so[["percent.mt"]])
summary(so[["cov1227"]])
summary(so[["nFeature_RNA"]])
VlnPlot(so, features = c("percent.mt", "cov1227", "nFeature_RNA"))

# filter seurat object
so                      <- subset(so, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 15 &
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
FeaturePlot(so, "af1227")

saveRDS(so, file = "~/Desktop/reznik/rrna-hotspots/hek_cells/data/hek_processed_seurat.rds")