rm(list=ls())
setwd("~/Desktop/reznik/rrna-hotspots/analysis_revisions/")

source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')

hotspots        <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t")
hotspots$label <- ifelse(hotspots$mutations_at_pos > 25 | hotspots$pos == 1227, hotspots$pos, NA)

panel_1b <- ggplot() +
  geom_segment(data = hotspots, aes(x = pos, xend = pos, y = 0, yend = mutations_at_pos), color = "#79706E", linewidth = 0.25) +
  geom_point(data = hotspots, aes(x = pos, y = mutations_at_pos, size = 0.5, alpha = 0.5, stroke = NA, colour = symbol), size = 2.5) + 
  scale_fill_identity() +
  theme_classic(base_size =  7,
                base_family = "ArialMT") + 
  scale_x_continuous(breaks = seq(min(hotspots$pos), max(filtered_gene_list$end), by = 1250), expand=(c(0,620))) +
  scale_y_continuous(expand = expansion(c(0,0.1)), n.breaks = 10) + 
  labs(x = "", 
       y = "Number of mutations at position") +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank()) +
  geom_text(data = hotspots, aes(x = pos, y = mutations_at_pos, label = label), 
            vjust = -1.5, size = 2, family = "ArialMT",fontface = "bold" , na.rm = TRUE) +
  scale_colour_manual(values = c("MT-RNR1" = "#59A14F", "MT-TV" =  "#86BCB6","MT-RNR2" = "#F1CE63", "MT-TL1" = "#86BCB6", "MT-ND1" = "#D7B5A6", "MT-TQ" = "#86BCB6",
                                 "MT-ND2" = "#C4C0D3","MT-TA" = "#86BCB6", "MT-CO1" = "#E15759", "MT-CO2" = "#D48B8C","MT-ATP8" = "#4A9894", 
                                 "MT-ATP8" ="#507AA8","MT-CO3"= "#D6B181", 
                                 "MT-ND3" = "#F28F28","MT-ND4L" = "#B69930",
                                "MT-ND4" =  "#8FCB7E", "MT-TL2" ="#4A9894", "MT-ND5" = "#9D7661", "MT-ND6" = "#BBB1AD", "MT-CYB" = "#FABFD2"  , "MT-TT" ="#4A9894" ))

ggsave(panel_1b, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig1b_hotspots.pdf", width = 4, height = 2)