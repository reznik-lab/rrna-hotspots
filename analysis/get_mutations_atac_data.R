## get variants 
gc()
setwd("~/Desktop/reznik/rrna-hotspots/analysis/")
source("./variant_calling.R")

## Adapted from Caleb Lareau ####

BErep1all        <- readRDS("~/Desktop/reznik/1227_atac/data/BErep1.rds")
BErep2all        <- readRDS("~/Desktop/reznik/1227_atac/data/BErep2.rds")
WTall            <- readRDS("~/Desktop/reznik/1227_atac/data/WT.rds")

BErep1MUTS       <- call_mutations_mgatk(BErep1all, stabilize_variance = FALSE)
BErep2MUTS       <- call_mutations_mgatk(BErep2all, stabilize_variance = FALSE) 
WTMUTS           <- call_mutations_mgatk(WTall, stabilize_variance = FALSE) 

BE1_all_muts     <-  data.frame(rowData(BErep1MUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)
BE2_all_muts     <-  data.frame(rowData(BErep2MUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)
WT_all_muts      <-  data.frame(rowData(WTMUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)

BE1_af           <- as.data.frame(assays(BErep1MUTS)[["allele_frequency"]])
BE1_cov          <- as.data.frame(assays(BErep1MUTS)[["coverage"]])

BE2_af           <- as.data.frame(assays(BErep2MUTS)[["allele_frequency"]])
BE2_cov          <- as.data.frame(assays(BErep2MUTS)[["coverage"]])

WT_af            <- as.data.frame(assays(WTMUTS)[["allele_frequency"]])
WT_cov           <- as.data.frame(assays(WTMUTS)[["coverage"]])

get_specific_var <- function(af, cov, mutation, filename = NA){
  position       <- substr(mutation, 1, nchar(mutation) -3)
  
  af             <- as.data.frame(t(af[mutation, ]))
  cov            <- as.data.frame(t(cov[mutation, ]))
  
  af$rows        <- rownames(af)
  cov$rows       <- rownames(cov)
  df             <- merge(af, cov, by = "rows")
  rownames(df)   <- df$cell
  colnames(df) <- c("cell", paste0("af",position), paste0("cov", position))
  
  if(!is.na(filename)){
    write.csv(df, file = filename, row.names = FALSE)
  }
  
  
  return(df)
}

BE1_1227 <- get_specific_var(BE1_af, BE1_cov, "1227G>A") 
BE1_9369 <- get_specific_var(BE1_af, BE1_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-9369.csv")
BE1_1230 <- get_specific_var(BE1_af, BE1_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1230.csv")
BE1_1233 <- get_specific_var(BE1_af, BE1_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1233.csv")

merged_BE1 <- merge(BE1_1227, BE1_9369, by = "cell")
merged_BE1 <- merge(merged_BE1, BE1_1230, by = "cell") 
merged_BE1 <- merge(merged_BE1, BE1_1233, by = "cell") 

merged_BE1 <- subset(merged_BE1, cov1227 > 10 & cov1230 > 10 & cov1233 > 10 & cov9369 > 10)
WT_BE1 <- merged_BE1 %>% subset(af1230 < 0.05 & af1233 < 0.05 & af9369 < 0.05 & af1227 > 0.01)
WT_cells <- WT_BE1$cell
hist(WT_BE1$af1227)
BE2_1227 <- get_specific_var(BE2_af, BE2_cov, "1227G>A")
BE2_9369 <- get_specific_var(BE2_af, BE2_cov, "9369T>C")

BE2_9369 <- get_specific_var(BE2_af, BE2_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-9369.csv")
BE2_1230 <- get_specific_var(BE2_af, BE2_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-1230.csv")
BE2_1233 <- get_specific_var(BE2_af, BE2_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-1233.csv")

WT_9369 <- get_specific_var(WT_af, WT_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/WT-9369.csv")
WT_1230 <- get_specific_var(WT_af, WT_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/WT-1230.csv")
WT_1233 <- get_specific_var(WT_af, WT_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/WT-1233.csv")

merged_BE2 <- merge(BE2_1227, BE2_9369, by = "cell")

ggplot(merged_BE1, aes(x = af1227, y = af9369)) + 
  geom_point(alpha = 0.5)

rm(BErep1all, BErep2all, WTall, BErep1MUTS, BErep2MUTS, WTMUTS)
gc()
