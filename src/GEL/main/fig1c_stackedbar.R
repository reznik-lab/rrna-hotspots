## AUTHOR: SONIA BOSCENCO

rm(list=ls())
source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')

gel_annotations        <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t")
gel_annotations$class  <- ifelse(grepl("^MT-T", gel_annotations$symbol), "tRNA",
                                 ifelse(grepl("^MT-R", gel_annotations$symbol), "rRNA", "Protein Coding"))

n_genes_pc             <- length(unique(gel_annotations[gel_annotations$class == "Protein Coding", ]$symbol))
n_genes_rRNA           <- length(unique(gel_annotations[gel_annotations$class == "rRNA", ]$symbol))
n_genes_tRNA           <- length(unique(gel_annotations[gel_annotations$class == "tRNA", ]$symbol))

table_annotations      <- as.data.frame(table(gel_annotations$class))
table_annotations$prop <- table_annotations$Freq / sum(table_annotations$Freq)
table_annotations$Var1 <- factor(table_annotations$Var1, levels = c("rRNA","tRNA","Protein Coding"))

use_colors             <- c("Protein Coding" = "#A0CBE8","tRNA" = "#86BCB6","rRNA" = "#499894")

panel_1c <- ggplot(table_annotations, aes(x = "", y = prop, fill = Var1)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.25, color = "white", linewidth = 0.1) + 
  theme_classic(base_family = "ArialMT",
                base_size = 7) +  
  scale_fill_manual(values = use_colors) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "", 
       y = 'Proportion of hotspots', 
       fill = "") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5, family = "ArialMT"),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(1, "pt"),
        axis.line.x = element_blank())

ggsave(panel_1c, file = '~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig1c_hotspotspropbarplot.pdf',width=1.75,height=2)