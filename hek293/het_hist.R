## AUTHOR: SONIA BOSCENCO
## create heteroplasmy violin plot

library(ggbeeswarm)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load heteroplasmy
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
het_df                 <- as.data.frame(rbind(fread("~/Desktop/reznik/rrna-hotspots/hek_cells/data/Experiment3_HEK293_1227G_A.csv")%>%
                                                 mutate(lib = "HEK293")))
                                        
                                    

het_df2                 <- as.data.frame(rbind(fread("~/Desktop/reznik/1227_atac//data/het/BErep1-1227.csv") %>%
                                                mutate(cell = gsub("-1", "-1", cell)) %>% mutate(lib = "143b Rep1"),
                                              fread("~/Desktop/reznik/1227_atac//data/het/BErep2-1227.csv") %>%
                                                mutate(cell = gsub("-1", "-2", cell)) %>% mutate(lib = "143b Rep2")))

het_df_lib1               <- het_df2[het_df2$lib == "143b Rep1", ]
het_df_lib1               <- het_df_lib1[het_df_lib1$cov >= 25, ]
het_df_lib1$af1227        <- as.numeric(het_df_lib1$af1227)

het_df_lib2               <- het_df2[het_df2$lib == "143b Rep2", ]
het_df_lib2               <- het_df_lib2[het_df_lib2$cov >= 25, ]
het_df_lib2$af1227        <- as.numeric(het_df_lib2$af1227)

het_df_hek                <- het_df[het_df$cov1227 >= 25, ]
het_df_hek$af1227         <- as.numeric(het_df_hek$af1227)

median_lib1               <- median(het_df_lib1$af1227)
median_lib2               <- median(het_df_lib2$af1227)
median_hek                <- median(het_df_hek$af1227)

mean(het_df_lib1$af1227)
mean(het_df_lib2$af1227)

perc_cells_less10_lib1    <- sum(het_df_lib1$af1227 < 0.1) / nrow(het_df_lib1)
perc_cells_less10_lib2    <- sum(het_df_lib2$af1227 < 0.1) / nrow(het_df_lib2)
perc_cells_less10_hek     <- sum(het_df_hek$af1227 < 0.1) / nrow(het_df_hek)

perc_cells_greater90_lib1 <- sum(het_df_lib1$af1227 > 0.80) / nrow(het_df_lib1)
perc_cells_greater90_lib2 <- sum(het_df_lib2$af1227 > 0.80) / nrow(het_df_lib2)
perc_cells_greater90_hek  <- sum(het_df_hek$af1227 > 0.80) / nrow(het_df_hek)

wt <- as.data.frame(fread("~/Desktop/reznik/1227_atac//data/het/WT-1227.csv") %>%
  mutate(cell = gsub("-1", "-3", cell)) %>% mutate(lib = "143b WT"))

het_df_all <- rbind(het_df, het_df2)
het_df_all <- subset(het_df_all, af1227 > 0)
het_df_all <- rbind(het_df_all, wt)
rownames(het_df)<- het_df$cell
hist                      <- ggplot(het_df, aes(x = af1227)) + 
  geom_histogram(bins = 50, linewidth = 0.01, fill = "black") + 
  theme_classic(base_size = 7, 
                base_family = "ArialMT") + 
  labs(x = "HEK m.1227G>A heteroplasmy", y = "Frequency") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"))

ggsave(hist, file = "~/Desktop/reznik/rrna-hotspots/hek_cells/figures/hist_heteroplasmy.pdf", width = 3, height = 2)


het_violin <- ggplot(het_df_all, aes(x = factor(lib, levels = c("HEK293", "143b Rep1", "143b Rep2", "143b WT")), y = af1227, color = lib)) + 
  geom_boxplot(linewidth = 0.1, colour = "black", outlier.shape = NA) + 
  geom_quasirandom(alpha = 0.1, stroke = NA, size = 1) + 
  theme_classic(base_size = 7, 
                base_family = "ArialMT") + 
  labs(x = "", y = "m.1227G>A heteroplasmy") + 
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt")) + 
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#499894", "#B07AA1"))

sub <- subset(het_df_all, lib %in% c("HEK293", "143b Rep1"))
het_violin <- ggplot(sub, aes(x = factor(lib, levels = c("HEK293", "143b Rep1")), y = af1227, color = lib)) + 
  geom_boxplot(linewidth = 0.1, colour = "black", outlier.shape = NA) + 
  geom_quasirandom(alpha = 0.1, stroke = NA, size = 1) + 
  theme_classic(base_size = 7, 
                base_family = "ArialMT") + 
  labs(x = "", y = "m.1227G>A heteroplasmy") + 
  theme(legend.position = "none",
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt")) + 
  scale_color_manual(values = c("#4E79A7", "#B07AA1"))

ggsave(het_violin, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig3/fig3a_heteroplasmyviolin.pdf", width = 2.25, height = 2.25)
