source('~/Desktop/reznik/rrna-hotspots/revisions/scripts/prerequisites.R')
library(ggplot2)
library(ggpubr)
library(scales)
library(dplyr)
gel_hotspots                   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
gel_mutations                  <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
annotations                    <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/mito_rRNA_annotations_NL.tsv", sep = "\t")
rnr_muts_annotated             <- gel_mutations %>%
                                  subset(Hugo_Symbol %in% c("MT-RNR1", "MT-RNR2")) %>% 
                                           group_by(Start_Position) %>%
                                           distinct(Start_Position, .keep_all = TRUE)

rnr_muts_annotated$annotation  <- annotations$RNA_base_type[match(rnr_muts_annotated$Start_Position,annotations$POS)]
rnr_muts_annotated$is_hotspot  <- ifelse(rnr_muts_annotated$Start_Position %in% gel_hotspots$pos,"Hotspot","Not Hotspot")

all_muts_melt                   <- melt((table(rnr_muts_annotated$annotation,rnr_muts_annotated$is_hotspot)))
colnames(all_muts_melt)         <- c("Annotation","Hotspot","Frequency")

all_muts_melt$Proportion        <- c(prop.table(all_muts_melt[all_muts_melt$Hotspot=="Hotspot",3]),prop.table(all_muts_melt[all_muts_melt$Hotspot=="Not Hotspot",3]))

use_cols                        <- c("#4E79A7", "#499894", "#B07AA1")
p.value                          <- fisher.test(table(rnr_muts_annotated$annotation, rnr_muts_annotated$is_hotspot))
panel_2e                        <- ggplot(all_muts_melt, aes(x=Hotspot, y=Proportion, fill=Annotation)) + 
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) + 
  theme_classic(base_size = 7,
                base_family = "ArialMT") + 
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        legend.key.size = unit(6, "pt")) +
  scale_fill_manual(values = use_cols) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) 


 #0.0005267
p.value$data
ggsave(panel_2e, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig2/fig2e_wc_barplot.pdf", width = 1.5, height = 1.5)

table_of_vals <- (table(rnr_muts_annotated$annotation, rnr_muts_annotated$is_hotspot))
wc_hotspot <- table_of_vals["WC", "Hotspot"]
non_wc_and_loop_hotspot <- table_of_vals["non-WC", "Hotspot"] + table_of_vals["loop-or-other", "Hotspot"]
wc_non_hotspot <- table_of_vals["WC", "Not Hotspot"]
non_wc_and_loop_non_hotspot <- table_of_vals["non-WC", "Not Hotspot"] + table_of_vals["loop-or-other", "Not Hotspot"]
OR_combined <- (wc_hotspot * non_wc_and_loop_non_hotspot) / (non_wc_and_loop_hotspot * wc_non_hotspot)
log2(OR_combined)
