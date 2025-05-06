rm(list=ls())
setwd("~/Desktop/reznik/rrna-hotspots/analysis_revisions/")

source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')

hotspots        <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/supp_tables/table_3_indel_hotspots.txt", sep = " ")


panel_2_supp <- ggplot() +
  geom_segment(data = hotspots, aes(x = start, xend = start, y = 0, yend = mutant), color = "#79706E", linewidth = 0.25) +
  geom_point(data = hotspots, aes(x = start, y = mutant, size = 0.5, alpha = 0.5, stroke = NA, colour = gene), size = 2.5) + 
  scale_fill_identity() +
  theme_classic(base_size =  7,
                base_family = "ArialMT") + 
  scale_x_continuous(breaks = seq(709, max(filtered_gene_list$end), by = 1250), expand=(c(0,620))) +
  scale_y_continuous(expand = expansion(c(0,0.1)), n.breaks = 10) + 
  labs(x = "", 
       y = "Number of mutations at position") +
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1)) +
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.line.x = element_blank()) +
  geom_text(data = hotspots, aes(x = start, y = mutant, label = homopolymer), 
            vjust = -1.5, size = 2, family = "ArialMT",fontface = "bold" , na.rm = TRUE) +
  scale_colour_manual(values = c("MT-ND1" = "#D7B5A6","MT-CO3"= "#D6B181", 
                                 "MT-ND3" = "#F28F28",
                                 "MT-ND4" =  "#8FCB7E", "MT-ND5" = "#9D7661"))

ggsave(panel_2_supp, file = "~/Desktop/reznik/rrna-hotspots/revisions/supplemental/fig2a_indel_hotspots.pdf", width = 4, height = 2)