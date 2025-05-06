be1ranks <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/results/be1_ranks.csv")
hekranks <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/results/hek_ranks.csv")

combined_ranks <- merge(be1ranks, hekranks, by = "gene")
combined_ranks$diff <- abs(combined_ranks$z.x - combined_ranks$z.y)
combined_ranks$colour <- ifelse(combined_ranks$gene  %in% c("AC093512.2","MT-CO1", "MT-CO3", "MT-ATP6", "MT-CYB", "RASL10A","AC011447.3", "ADIRF", "MT-CO2", "MT-ND4", "MT-ND5", "RPL35", "CAMK2N2", "MT-ND3"), TRUE, FALSE)
p <-ggplot(combined_ranks, aes(x = z.x, y = z.y, colour = colour)) + 
  geom_point(size = 1, alpha = 0.5, stroke = NA) + 
  labs(x = "143b gene ranks", y = 'HEK293 gene ranks') +
  theme_classic(base_size = 7, 
                base_family = 'ArialMT') + 
  theme(
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(1, "pt"),
        legend.position = "none") +
  scale_colour_manual(values = c("black", "#E15759")) + 
  geom_text_repel(aes(label = ifelse(gene %in% c("AC093512.2","MT-CO1", "MT-ATP6", "MT-CYB", "RASL10A","AC011447.3", "ADIRF", "MT-CO2", "MT-ND4", "MT-ND5", "RPL35", "CAMK2N2"), gene, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) 

ggsave(p, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/be1_v_hek_ranks.pdf", width = 2, height = 2.5, useDingbats = FALSE, dpi = 300)