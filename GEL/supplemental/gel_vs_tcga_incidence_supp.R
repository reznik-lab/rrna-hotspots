## AUTHOR: SONIA BOSCENCO

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparing incidence of mutations in PCAWG vs. GEL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R")
gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
tcga_mutations                 <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data//data_mutations_tcga.txt", sep = "\t", header = TRUE) %>% 
                                  subset(Variant_Type == "SNP")
tcga_hotspots                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data//hotspots_snv_precalculated.txt", sep = "\t", header = TRUE) %>%
                                  subset(q.value < 0.05)

# Have to figure out what to do with the masked counts, but none of these are hotpots anyways so will just
# set them all to 1 instead 
# Collapse gel mutations by positions 

gel_mutations$Freq         <- ifelse(gel_mutations$Freq == "< 5", 1, gel_mutations$Freq)
gel_mutations$Freq         <- as.numeric(gel_mutations$Freq)
gel_muts                   <- gel_mutations %>%
                              select(Start_Position, Freq) %>%
                              group_by(Start_Position) %>%
                              summarise(n = sum(Freq))

tcga_muts                  <- as.data.frame(table(tcga_mutations$Start_Position))

# Annotate hotspots
gel_muts$is_hotspot      <- ifelse(gel_muts$Start_Position %in% gel_hotspots$pos, TRUE, FALSE)
tcga_muts$is_hotspot     <- ifelse(tcga_muts$Var1 %in% tcga_hotspots$Start_Position, TRUE, FALSE)

# Merging
merged_tcga_gel          <- merge(gel_muts, tcga_muts, by.x = "Start_Position", by.y = "Var1", all = TRUE)
merged_tcga_gel[is.na(merged_tcga_gel)] <- 0

colnames(merged_tcga_gel)<- c("Start_Position", "GEL", "Hot_Gel", "TCGA", "Hot_TCGA")

merged_tcga_gel$col      <- ifelse(merged_tcga_gel$Hot_Gel == TRUE & merged_tcga_gel$Hot_TCGA == TRUE, "Hotspot GEL+TCGA",
                                    ifelse(merged_tcga_gel$Hot_Gel == FALSE & merged_tcga_gel$Hot_TCGA == TRUE, "Hotspot TCGA",
                                           ifelse(merged_tcga_gel$Hot_Gel == TRUE & merged_tcga_gel$Hot_TCGA == FALSE, "Hotspot GEL",  "Not Hotspot")))

# Ensure neither dots are behind the rest of the other dots
merged_tcga_gel           <- merged_tcga_gel[order(merged_tcga_gel$col != "Not Hotspot"), ]
color_order               <- c("Not Hotspot", "Hotspot GEL", "Hotspot TCGA", "Hotspot GEL+TCGA")
merged_tcga_gel$col       <- factor(merged_tcga_gel$col, levels = color_order)

supp2c                    <- ggplot(merged_tcga_gel, aes(x=GEL, y=TCGA, color=col)) + 
  geom_point(alpha = 0.9) +
  theme_classic() +
  stat_cor() +
  theme_std() + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks.length = unit(1, "pt")) + 
  scale_y_continuous(expand = c(0,1)) + 
  scale_x_continuous(expand = c(0,1)) +
  theme(legend.title = element_blank()) +
  xlab("Number samples mutated in GEL") + 
  ylab("Number samples mutated in PCAWG") +
  scale_color_manual(values = c("#E15759", "#4E79A7","#59A14F", "#BAB0AC"),
                     breaks = c("Hotspot GEL+TCGA", "Hotspot GEL", "Hotspot TCGA", "Not Hotspot"))

ggsave(supp2c, file = '~/Desktop/reznik/rrna-hotspots/revisions/supplemental/gel_vs_tcga.pdf',width=3.5,height=3)