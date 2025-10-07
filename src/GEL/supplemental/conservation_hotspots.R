source('~/Desktop/reznik/rrna-hotspots/revisions/scripts/prerequisites.R')
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
hotspots       <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t", header = TRUE)
maf            <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/maf_snp_5.txt", sep = "\t", header = TRUE)
conservation   <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/conservation_rate.tsv", sep = "\t", header = TRUE)

maf$ishot      <- ifelse(maf$Start_Position %in% hotspots$pos, TRUE, FALSE)
maf$conservation <- conservation$rate[match(maf$Start_Position, conservation$pos)]
maf <- maf %>% distinct(Start_Position, .keep_all = TRUE)
rnr_maf      <- maf %>% subset(grepl("MT-R", Hugo_Symbol))

pc_maf        <- maf %>% subset(!grepl("MT-R", Hugo_Symbol) & !grepl("MT-T", Hugo_Symbol))
p <- ggplot(rnr_maf, aes(x = as.factor(ishot), y = conservation, fill = ishot)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(alpha = 0.1, stroke = 0, orientation = "x", bandwidth = 1) +
  theme_bw(base_family = "ArialMT", 
                base_size = 7) + 
  labs(x = "", y = "Conservation Rate [%]") + 
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") + 
  scale_x_discrete(labels = c("FALSE" = "Non-Hotspot", "TRUE" = "Hotspot"))+ 
  stat_compare_means(family = "ArialMT", size = 2,
                     tip.length = 0, label = "p.format", label.y = 105, label.x =1.5,
                     vjust =0, bracket.size = 0.1) + 
 
  scale_fill_manual(values = c("#A6CEE3", "#b2df8a"))

ggsave(p, file = "~/Desktop/reznik/rrna-hotspots/revisions/results/rrna_hots_conservationboxplot.pdf", width = 3, height = 3)
ggplot(pc_maf, aes(x = ishot, y = conservation, fill = ishot)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.1, stroke = 0) + 
  stat_compare_means(comparisons = list(c(FALSE, TRUE)), na.rm = TRUE, show.legend = FALSE, size = 1, tip.length = 0, position = position_nudge(y = 2)) + 
  theme_minimal(base_family = "ArialMT", 
           base_size = 7) + 
  labs(x = "", y = "Conservation Rate [%]") + 
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels = c("FALSE" = "Non-Hotspot", "TRUE" = "Hotspot")) + 
  scale_fill_manual(values = c("#A6CEE3", "#b2df8a"))