## get variants 
gc()
setwd("~/Desktop/reznik/rrna-hotspots/hek_cells/analysis")
source("./variant_calling.R")

## Adapted from Caleb Lareau ####
## call mutations in hek cells 


HEK              <- readRDS("~/Desktop/reznik/rrna-hotspots/hek_cells/data/Experiment3_HEK293.rds")
muts             <- call_mutations_mgatk(HEK, stabilize_variance = FALSE)
all_muts         <- data.frame(rowData(muts)) %>% filter(n_cells_conf_detected > 10) %>% filter(mean > 0.05)

af               <- as.data.frame(assays(muts)[["allele_frequency"]])
cov              <- as.data.frame(assays(muts)[["coverage"]])

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

BE1_1227 <- get_specific_var(af, cov, "1227G>A") 
BE1_9369 <- get_specific_var(af, cov, "9369T>C", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-9369.csv")
BE1_1230 <- get_specific_var(af, cov, "1230C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1230.csv")
BE1_1233 <- get_specific_var(af, cov, "1233C>T", filename = "~/Desktop/reznik/1227_atac//data/het/BErep1-1233.csv")
