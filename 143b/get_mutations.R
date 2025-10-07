## AUTHOR: SONIA BOSCENCO
## get all variants including off target

gc()
source("~/Desktop/reznik/1227_atac/analysis/variant_calling.R")
source("~/Desktop/reznik/1227_atac/prerequisites.R")

# load data 
BErep1all        <- readRDS("~/Desktop/reznik/1227_atac/data/BErep1.rds")
BErep2all        <- readRDS("~/Desktop/reznik/1227_atac/data/BErep2.rds")
WTall            <- readRDS("~/Desktop/reznik/1227_atac/data/WT.rds")

# call mutations using mgatk
BErep1MUTS       <- call_mutations_mgatk(BErep1all, stabilize_variance = FALSE)
BErep2MUTS       <- call_mutations_mgatk(BErep2all, stabilize_variance = FALSE) 
WTMUTS           <- call_mutations_mgatk(WTall, stabilize_variance = FALSE) 

BE1_all_muts     <-  data.frame(rowData(BErep1MUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)
BE2_all_muts     <-  data.frame(rowData(BErep2MUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)
WT_all_muts      <-  data.frame(rowData(WTMUTS)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)

# save info 
write.csv(BE1_all_muts, file = "~/Desktop/reznik/rrna-hotspots/data/processed/SI_figures/BE1_all_muts.csv", row.names = FALSE)
write.csv(BE2_all_muts, file = "~/Desktop/reznik/rrna-hotspots/data/processed/SI_figures/BE2_all_muts.csv", row.names = FALSE)
write.csv(WT_all_muts, file = "~/Desktop/reznik/rrna-hotspots/data/processed/SI_figures/WT_all_muts.csv", row.names = FALSE)

BE1_all_muts[order(BE1_all_muts$mean),]$position
ggplot(BE1_all_muts, aes(x = factor(position, levels = BE1_all_muts[order(BE1_all_muts$mean),]$position), y = mean)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_classic() + 
  labs(x = "") + 
  scale_y_continuous(expand = c(0,0))

# recover allele frequency and coverage
BE1_af           <- as.data.frame(as.matrix(assays(BErep1MUTS)[["allele_frequency"]]))
BE1_cov          <- as.data.frame(as.matrix(assays(BErep1MUTS)[["coverage"]]))

BE2_af           <- as.data.frame(as.matrix(assays(BErep2MUTS)[["allele_frequency"]]))
BE2_cov          <- as.data.frame(as.matrix(assays(BErep2MUTS)[["coverage"]]))

WT_af            <- as.data.frame(as.matrix(assays(WTMUTS)[["allele_frequency"]]))
WT_cov           <- as.data.frame(as.matrix(assays(WTMUTS)[["coverage"]]))

# function to get info for specific variants
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

# choose vars of interest 
BE1_1227 <- get_specific_var(BE1_af, BE1_cov, "1227G>A") 
BE1_9369 <- get_specific_var(BE1_af, BE1_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-9369.csv")
BE1_1230 <- get_specific_var(BE1_af, BE1_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1230.csv")
BE1_1233 <- get_specific_var(BE1_af, BE1_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1233.csv")

BE2_1227 <- get_specific_var(BE2_af, BE2_cov, "1227G>A")
BE2_9369 <- get_specific_var(BE2_af, BE2_cov, "9369T>C")

BE2_9369 <- get_specific_var(BE2_af, BE2_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-9369.csv")
BE2_1230 <- get_specific_var(BE2_af, BE2_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-1230.csv")
BE2_1233 <- get_specific_var(BE2_af, BE2_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep2-1233.csv")

WT_9369 <- get_specific_var(WT_af, WT_cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/WT-9369.csv")
WT_1230 <- get_specific_var(WT_af, WT_cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/WT-1230.csv")
WT_1233 <- get_specific_var(WT_af, WT_cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/WT-1233.csv")


