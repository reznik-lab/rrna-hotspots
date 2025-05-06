
maf                 <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t")
pcawg_mutations     <- read.csv("~/Desktop/reznik/rrna-hotspots/data/processed/fig1_2_tables/data_mutations_pcawg.txt", sep = "\t", header = TRUE) %>% 
                       subset(Variant_Type == "SNP")
tcga_mutations      <- read.csv("~/Desktop/reznik/rrna-hotspots/data/processed/fig1_2_tables/data_mutations_tcga.txt", sep = "\t", header = TRUE) %>% 
                       subset(Variant_Type == "SNP")

table_pcawg        <- as.data.frame(table(pcawg_mutations$ShortVariantID))
table_tcga         <- as.data.frame(table(tcga_mutations$ShortVariantID))

maf$ShortVariantID <- str_replace_all(string = maf$ShortVariantID, pattern = ":", replacement = "")
maf$pcawg_freq      <- table_pcawg$Freq[match(maf$ShortVariantID, table_pcawg$Var1)]
maf$tcga_freq      <- table_tcga$Freq[match(maf$ShortVariantID, table_tcga$Var1)]

write.table(maf, file = "~/Desktop/reznik/rrna-hotspots/revisions/supp_tables/gel_collapsed_maf.txt", sep = "\t", row.names = FALSE)

