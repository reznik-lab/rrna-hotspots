source('~/Desktop/reznik/rrna-hotspots/revisions/scripts/prerequisites.R')

gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
helix_muts          <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/helix_mutations.csv", sep = ",", header = TRUE, row.names=1)  %>% 
  subset(variant_type == "SNP") 

rnr_positions                  <- data.frame(c(648:1601, 1671:3229))
colnames(rnr_positions)        <- c('position')

helix_collapsed                <- helix_muts %>%
  subset(gene %in% c("MT-RNR1", "MT-RNR2")) %>% 
  group_by(Start_Position) %>% 
  summarize(total_counts = sum(total_count), total_homs = sum(counts_hom), total_hets = sum(counts_het))

helixs_conserved               <- merge(rnr_positions, helix_collapsed,by.x='position',by.y='Start_Position',all.x = T)
helixs_conserved[(is.na(helixs_conserved))] <- 0

# subsetting for only RNR genes
rnr_hotpots                    <- gel_hotspots %>% 
  subset(symbol %in% c("MT-RNR1", "MT-RNR2"))

rnr_muts                       <- gel_mutations %>% 
  subset(Hugo_Symbol %in% c("MT-RNR1", "MT-RNR2"))

# classify all mutations as hotpots or not
helixs_conserved$is_hotspot   <- ifelse(helixs_conserved$position %in% rnr_hotpots$pos, "Hotspot", "Not Hotspot")

# can remove positions that have more than one mutation
rnr_muts                       <- rnr_muts %>%
  group_by(Start_Position) %>%
  distinct(Start_Position, .keep_all = TRUE)

helixs_conserved               <- merge(helixs_conserved ,rnr_muts, by.x='position', by.y='Start_Position', all.y = T)

helixs_conserved$is_conserved  <- ifelse(helixs_conserved$total_homs < 2, "Constrained", "Not Constrained")

table_of_constrained_hots      <- as.data.frame(table(helixs_conserved$is_conserved, helixs_conserved$is_hotspot))

n_hots                         <- sum(table_of_constrained_hots[table_of_constrained_hots$Var2 == "Hotspot", ]$Freq)
n_not_hots                     <- sum(table_of_constrained_hots[table_of_constrained_hots$Var2 == "Not Hotspot", ]$Freq)
table_of_constrained_hots$prop <- table_of_constrained_hots$Freq / c(n_hots, n_hots, n_not_hots,n_not_hots)
table_of_constrained_hots$Var1 <- factor(table_of_constrained_hots$Var1, levels = c("Not Constrained", "Constrained"))

use_cols                       <- c("#BAB0AC", "#79706E")

panel_2h                       <- ggplot(table_of_constrained_hots, aes(x = Var2, y = prop, fill = Var1)) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = "white") + 
  labs(x = "Hotspot",
       y = "Proportion",
       fill = "") + 
  scale_fill_manual(values = use_cols) +
  theme_classic(base_size = 7, 
                base_family = "ArialMT") + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        legend.key.size = unit(6, "pt")) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))

ggsave(panel_2h, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig2/fig2h_constrained_barplot.pdf", width = 1.5, height = 1.5)


p.value     <-  fisher.test(matrix(table_of_constrained_hots$Freq, ncol = 2)) #0.04908
p.value




helixs_conserved$het_conserved <- ifelse(helixs_conserved$total_hets < 2, "Constrained", "Not Constrained")

table_of_constrained_hets      <- as.data.frame(table(helixs_conserved$het_conserved, helixs_conserved$is_hotspot))

n_hots                         <- sum(table_of_constrained_hets[table_of_constrained_hets$Var2 == "Hotspot", ]$Freq)
n_not_hots                     <- sum(table_of_constrained_hets[table_of_constrained_hets$Var2 == "Not Hotspot", ]$Freq)
table_of_constrained_hets$prop <- table_of_constrained_hets$Freq / c(n_hots, n_hots, n_not_hots,n_not_hots)
table_of_constrained_hets$Var1 <- factor(table_of_constrained_hets$Var1, levels = c("Not Constrained", "Constrained"))

use_cols                       <- c("#499894", "#B07AA1")

p <- ggplot(table_of_constrained_hets, aes(x = Var2, y = prop, fill = Var1)) + 
  geom_bar(stat = "identity", linewidth = 0.1, colour = "white") + 
  labs(x = "Hotspot",
       y = "Proportion",
       fill = "") + 
  scale_fill_manual(values = use_cols) +
  theme_classic(base_size = 7, 
                base_family = "ArialMT") + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        legend.key.size = unit(6, "pt")) + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))

fisher.test(matrix(table_of_constrained_hets$Freq, ncol = 2))

ggsave(p, file = "~/Desktop/reznik/rrna-hotspots/revisions/supplemental/heteroplasmy_constrained_barplot.pdf", width = 1.5, height = 1.5)