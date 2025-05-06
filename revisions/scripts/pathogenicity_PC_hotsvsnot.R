# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare pathogenecity scores of PC hotspots to non
# hotspots 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')
gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
pathogenic_pc       <- read.csv("~/Desktop/reznik/rrna-hotspots/data/processed/fig1_2_tables/MitImpact_db_3.1.2.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
helix_muts          <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/helix_mutations.csv", sep = ",", header = TRUE, row.names=1)  %>% 
  subset(variant_type == "SNP") 
# get only protein coding genes
pc_genelist                    <- gene_list[!grepl("^R|^T|CR", gene_list$gene), ]

pc_muts                        <- gel_mutations[!grepl("^MT-R|^MT-T|CR|ControlRegion|Translation_Start_Site", gel_mutations$Hugo_Symbol), ]
pc_muts                        <- pc_muts[!grepl("Translation_Start_Site", pc_muts$Variant_Classification), ]

# can remove positions that have more than one mutation
pc_muts                        <- pc_muts %>%
  group_by(Start_Position) %>%
  distinct(Start_Position, .keep_all = TRUE)

pc_muts$is_hotspot <- ifelse(pc_muts$Start_Position %in% pc_hotspots$pos, "Hotspot", "Not Hotspot")
pc_muts$Alt <- str_sub(pc_muts$ShortVariantID, -1, -1)
pc_pathogenics <- merge(pc_muts, pathogenic_pc, by.x = c("Start_Position", "Alt"), by.y = c("Start", "Alt"), all.x = TRUE)

pc_pathogenics$APOGEE2 <- ifelse(pc_pathogenics$APOGEE2 %in% c("VUS+", "VUS-"), "VUS", pc_pathogenics$APOGEE2)
table_scores <- as.data.frame(table(pc_pathogenics$APOGEE2, pc_pathogenics$is_hotspot))

n_pc_muts <- sum(table_scores[table_scores$Var2 == "Not Hotspot",]$Freq)
n_pc_hots <- sum(table_scores[table_scores$Var2 == "Hotspot",]$Freq)

table_scores$n <- ifelse(table_scores$Var2 == "Hotspot", n_pc_hots, n_pc_muts)
table_scores$prop <- table_scores$Freq / table_scores$n

table_scores$Var1 <- factor(table_scores$Var1, levels = rev(c("VUS","Benign", "Likely-benign", "Likely-pathogenic", "Pathogenic")))

use_cols                       <- c("#bebada", "#8dd3c7", "#d9d9d9", "#d9d9d9", "#d9d9d9")
patho <- ggplot(table_scores, aes(x = Var2, y = prop, fill = Var1)) + 
  geom_bar(stat = "identity", size = 0.2, colour = "white") +  
  scale_fill_manual(values = use_cols) +
  theme_std() + 
  labs(x = "",
       y = "Proportion",
       fill = "APOGEE2 Classification")

if(get_pdfs==T) ggsave('../figures_v2/panel_2l.pdf',width=2.5,height=2)
p.value                          <- fisher.test(table(pc_pathogenics$APOGEE2, pc_pathogenics$is_hotspot)) #0.004478
table <- table(pc_pathogenics$APOGEE2, pc_pathogenics$is_hotspot)
# just trying to compare the idfference between patho / liekly patho categories and everything else 
combined <- rbind(
  "Likely/Pathogenic" = colSums(table[c("Likely-pathogenic", "Pathogenic"), ]),
  "Other" = colSums(table[c("Benign", "Likely-benign", "VUS", "VUS-", "VUS+"), ])
)
fisher.test(matrix(combined, nrow=2))