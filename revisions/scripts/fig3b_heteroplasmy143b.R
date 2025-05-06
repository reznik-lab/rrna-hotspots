so <- readRDS("~/Desktop/reznik/rrna-hotspots/revisions/data/seurat_143b_processed.rds")
gex_mat                    <- Read10X_h5("~/Desktop/reznik/1227_atac//data/aggr/filtered_feature_bc_matrix.h5")$`Gene Expression`
all.genes            <- rownames(gex_mat)

correlate_het            <- function(gene_score_mat){
  ranks_list <- list()
  for (gene in colnames(gene_score_mat)) {
    if (gene == "het") next 
    
    nonzero_cells <- gene_score_mat[, gene] > 0
    
    #only test genes with more than 50 cells express it
    if (sum(nonzero_cells) > 50) {  
      cor_value <- cor(
        gene_score_mat[nonzero_cells, gene], 
        gene_score_mat[nonzero_cells, "het"], 
        method = "spearman", 
        use = "pairwise.complete.obs"
      )
      
      if (!is.na(cor_value)) {
        ranks_list[[gene]] <- cor_value
      }
    }
  }
  
  ranks <- sort(unlist(ranks_list), decreasing = TRUE)
  return(ranks)
}

gene_expression_mat      <- as.data.frame(t(so@assays$RNA$data))
het                      <- as.data.frame(so$af1227)
cov                      <- as.data.frame(so$cov1227)
cluster                  <- as.data.frame(so$seurat_clusters)
lib                      <- as.data.frame(so$lib)

colnames(het)            <- c("het")
colnames(cluster)        <- c("cluster")
colnames(cov)            <- c("cov")
colnames(lib)            <- c("lib")

gene_expression_mat$rows <- rownames(gene_expression_mat)
het$rows                 <- rownames(het)
cluster$rows             <- rownames(cluster)
cov$rows                 <- rownames(cov)
lib$rows                 <- rownames(lib)

master_df                <- merge(gene_expression_mat, het, by = "rows")
master_df                <- merge(master_df, cluster, by = "rows")
master_df                <- merge(master_df, cov, by = "rows")
master_df                <- merge(master_df, lib, by = "rows")


BE1                      <- subset(master_df, lib == "BErep1")
BE2                      <- subset(master_df, lib == "BErep2")

BE1_master_1227         <- BE1 %>% subset(het > 0.05)
BE2_master_1227         <- BE2 %>% subset(het > 0.05)

BE1_ranks               <- correlate_het(select(BE1_master_1227, -c("cluster", "rows", "cov", "lib")))
BE2_ranks               <- correlate_het(select(BE2_master_1227, -c("cluster", "rows", "cov", "lib")))

### mt CO1 correlation plot 

p1 <- ggplot(subset(BE1_master_1227, `MT-CO1` > 0), aes(x = het, y = `MT-CO1`)) + 
  geom_point(alpha = 0.5, colour = "#4E79A7", stroke = NA, size = 1) + 
  labs(x = "Heteroplasmy", y = "MT-CO1 Expression") + 
  scale_y_continuous(expand = c(0,0)) + 
  #geom_smooth(method = "lm", color = "black", linewidth = 0.2) + 
  stat_cor(family = "ArialMT",size=2, label.x.npc = "left", label.y.npc = "bottom") + 
  scale_x_continuous(limits = c(0,1)) + 
  
  theme_classic(base_size = 7, 
                base_family = 'ArialMT') +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("143b Rep 1") 

p2 <- ggplot(subset(BE2_master_1227, `MT-CO1` > 0), aes(x = het, y = `MT-CO1`)) + 
  geom_point(alpha = 0.5, colour = "#59A14F", stroke = NA, size = 1) + 
  labs(x = "Heteroplasmy", y = "MT-CO1 Expression") + 
  scale_y_continuous(expand = c(0,0)) + 
  #geom_smooth(method = "lm", color = "black", linewidth = 0.2) + 
  stat_cor(family = "ArialMT",size=2, label.x.npc = "left", label.y.npc = "bottom") + 
  scale_x_continuous(limits = c(0,1)) + 
  
  theme_classic(base_size = 7, 
                base_family = 'ArialMT') +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("143b Rep 2") 

ggsave(p1, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig3/fig3b_heteroplasmy_143b1.pdf", width = 2, height = 1.25)
ggsave(p2, file = "~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/fig3/fig3b_heteroplasmy_143b2.pdf", width = 2, height = 2)

