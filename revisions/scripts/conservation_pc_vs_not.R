# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Barplot comparing proportion of conserved vs. 
# non conserved positions in hotspots vs. not hotspots
# in PC genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')
gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
helix_muts          <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/helix_mutations.csv", sep = ",", header = TRUE, row.names=1)  %>% 
  subset(variant_type == "SNP") 
# get only protein coding genes
pc_genelist                    <- gene_list[!grepl("^R|^T|CR", gene_list$gene), ]

# create a list of all the positions in sequential order 
sequences_list                 <- apply(pc_genelist, 1, function(x) seq(from = x['start'], to = x['end']))
all_pos_pc                     <- unlist(sequences_list)
all_pos_pc_df                  <- data.frame(position = all_pos_pc)

# collapse helix with just cp

helix_pc                       <- helix_muts[!grepl("^MT-R|^MT-T|CR", helix_muts$gene), ]
helix_collapsed_pc             <- helix_pc %>%
  group_by(Start_Position) %>% 
  summarize(total_counts = sum(total_count), total_homs = sum(counts_hom), total_hets = sum(counts_het))

helixs_conserved_pc            <- merge(all_pos_pc_df, helix_collapsed_pc,by.x='position',by.y='Start_Position',all.x = T)
helixs_conserved_pc[(is.na(helixs_conserved_pc))] <- 0

pc_hotspots                    <- gel_hotspots[!grepl("^MT-R|^MT-T|CR|ControlRegion|Translation_Start_Site", gel_hotspots$symbol), ]
pc_muts                        <- gel_mutations[!grepl("^MT-R|^MT-T|CR|ControlRegion|Translation_Start_Site", gel_mutations$Hugo_Symbol), ]
pc_muts                        <- pc_muts[!grepl("Translation_Start_Site", pc_muts$Variant_Classification), ]

# classify all mutations as hotpots or not
helixs_conserved_pc$is_hotspot <- ifelse(helixs_conserved_pc$position %in% pc_hotspots$pos, "Hotspot", "Not Hotspot")

# can remove positions that have more than one mutation
pc_muts                        <- pc_muts %>%
  group_by(Start_Position) %>%
  distinct(Start_Position, .keep_all = TRUE)

# merge germline conservation with all mutations 
helixs_conserved_pc             <- merge(helixs_conserved_pc ,pc_muts, by.x='position', by.y='Start_Position', all.y = T)
helixs_conserved_pc$is_conserved<- ifelse(helixs_conserved_pc$total_homs < 2, "Constrained", "Not Constrained")

table_of_constrained_hots      <- as.data.frame(table(helixs_conserved_pc$is_conserved, helixs_conserved_pc$is_hotspot))

n_hots                         <- sum(table_of_constrained_hots[table_of_constrained_hots$Var2 == "Hotspot", ]$Freq)
n_not_hots                     <- sum(table_of_constrained_hots[table_of_constrained_hots$Var2 == "Not Hotspot", ]$Freq)
table_of_constrained_hots$prop <- table_of_constrained_hots$Freq / c(n_hots, n_hots, n_not_hots,n_not_hots)
table_of_constrained_hots$Var1 <- factor(table_of_constrained_hots$Var1, levels = c("Not Constrained", "Constrained"))

use_cols                       <- c("#ccebc5", "#bebada")

panel_supp                     <- ggplot(table_of_constrained_hots, aes(x = Var2, y = prop, fill = Var1)) + 
  geom_bar(stat = "identity") + 
  labs(x = "Hotspot",
       y = "Proportion",
       fill = "") + 
  scale_fill_manual(values = use_cols) +
  theme_std() + 
  theme(legend.position = "right",
        axis.title.x = element_blank())

if(get_pdfs==T) ggsave('../figures_v2/panel_2k.pdf',width=2.5,height=2)
p.value     <-  fisher.test(matrix(table_of_constrained_hots$Freq, ncol = 2)) #0.004554
p.value