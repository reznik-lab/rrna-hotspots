## AUTHOR: SONIA BOSCENCO


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparing incidence of mutations in PCAWG vs. GEL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R")
gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
pcawg_mutations                <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data//data_mutations_pcawg.txt", sep = "\t", header = TRUE) %>% 
                                  subset(Variant_Type == "SNP")
pcawg_hotspots                 <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data//hotspots_snv_pcawg.txt", sep = "\t", header = TRUE) %>%
                                  subset(q.value < 0.05)

pcawg_mutations <- pcawg_mutations %>% subset(Hugo_Symbol != "Control_Region")
# Have to figure out what to do with the masked counts, but none of these are hotpots anyways so will just
# set them all to 1 instead 
# Collapse gel mutations by positions 

gel_mutations$Freq         <- ifelse(gel_mutations$Freq == "< 5", 1, gel_mutations$Freq)
gel_mutations$Freq         <- as.numeric(gel_mutations$Freq)
gel_muts                   <- gel_mutations %>%
                              select(Start_Position, Freq) %>%
                              group_by(Start_Position) %>%
                              summarise(n = sum(Freq))

pcawg_muts               <- as.data.frame(table(pcawg_mutations$Start_Position))

# Annotate hotspots
gel_muts$is_hotspot      <- ifelse(gel_muts$Start_Position %in% gel_hotspots$pos, TRUE, FALSE)
pcawg_muts$is_hotspot    <- ifelse(pcawg_muts$Var1 %in% pcawg_hotspots$Start_Position, TRUE, FALSE)

# Merging
merged_pcawg_gel         <- merge(gel_muts, pcawg_muts, by.x = "Start_Position", by.y = "Var1", all = TRUE)
merged_pcawg_gel[is.na(merged_pcawg_gel)] <- 0

colnames(merged_pcawg_gel)<- c("Start_Position", "GEL", "Hot_Gel", "PCAWG", "Hot_PCAWG")

merged_pcawg_gel$col      <- ifelse(merged_pcawg_gel$Hot_Gel == TRUE & merged_pcawg_gel$Hot_PCAWG == TRUE, "Hotspot GEL+PCAWG",
                                    ifelse(merged_pcawg_gel$Hot_Gel == FALSE & merged_pcawg_gel$Hot_PCAWG == TRUE, "Hotspot PCAWG",
                                           ifelse(merged_pcawg_gel$Hot_Gel == TRUE & merged_pcawg_gel$Hot_PCAWG == FALSE, "Hotspot GEL",  "Not Hotspot")))

# Ensure neither dots are behind the rest of the other dots
merged_pcawg_gel          <- merged_pcawg_gel[order(merged_pcawg_gel$col != "Not Hotspot"), ]
color_order               <- c("Not Hotspot", "Hotspot GEL", "Hotspot PCAWG", "Hotspot GEL+PCAWG")
merged_pcawg_gel$col      <- factor(merged_pcawg_gel$col, levels = color_order)

supp2b                    <- ggplot(merged_pcawg_gel, aes(x=GEL, y=PCAWG, color=col)) + 
  #stat_cor() +
  geom_point(alpha = 0.9, stroke = NA) +
  theme_classic() +
  theme_std() + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks.length = unit(1, "pt")) + 
  scale_y_continuous(expand = c(0,1)) + 
  scale_x_continuous(expand = c(0,1)) +
  theme(legend.title = element_blank()) +
  xlab("Number samples mutated in GEL") + 
  ylab("Number samples mutated in PCAWG") +
  scale_color_manual(values = c("#E15759", "#4E79A7","#59A14F", "#BAB0AC"),
                     breaks = c("Hotspot GEL+PCAWG", "Hotspot GEL", "Hotspot PCAWG", "Not Hotspot"))

ggsave(supp2b, file = '~/Desktop/reznik/rrna-hotspots/revisions/supplemental/gel_vs_pcawg.pdf',width=3.5,height=3)