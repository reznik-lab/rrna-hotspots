## AUTHOR:: SONIA BOSCENCO

rm(list = ls())
library(Seurat)
so <- readRDS("~/Desktop/reznik/rrna-hotspots/revisions/data/seurat_143b_processed.rds")


gene_expression_mat      <- as.data.frame(t(so@assays$RNA$data))
het                      <- as.data.frame(so$af1227)
cov                      <- as.data.frame(so$cov1227)
cluster                  <- as.data.frame(so$seurat_clusters)
lib                      <- as.data.frame(so$lib)

colnames(het)            <- c("het")
colnames(cluster)        <- c("cluster")
colnames(lib)            <- c("lib")
colnames(cov)            <- c("cov")

gene_expression_mat$rows <- rownames(gene_expression_mat)
het$rows                 <- rownames(het)
cluster$rows             <- rownames(cluster)
lib$rows                 <- rownames(lib)
cov$rows                 <- rownames(cov)


master_df                <- merge(gene_expression_mat, het, by = "rows")
master_df                <- merge(master_df, cluster, by = "rows")
master_df                <- merge(master_df, lib, by = "rows")
master_df                <- merge(master_df, cov, by = "rows")

master_df                <- subset(master_df, het > 0.05)
BE1                      <- subset(master_df, lib == "BErep1")
BE2                      <- subset(master_df, lib == "BErep2")
WT                       <- subset(master_df, lib == "WT")


correlate_het_by_bins <- function(gene_score_mat) {
  gene_names <- colnames(gene_score_mat)[colnames(gene_score_mat) != "het"]
  results <- data.frame()
  
  for (bin_value in seq(0, 1.1, by = 0.1)) {
    above_bin <- gene_score_mat[gene_score_mat$het >= bin_value, ]
    below_bin <- gene_score_mat[gene_score_mat$het < bin_value, ]
    
    if (nrow(above_bin) > 1 & nrow(below_bin) > 1) {
      cor_above <- sapply(gene_names, function(gene) {
        nonzero_cells <- above_bin[[gene]] > 0
        if (sum(nonzero_cells, na.rm = TRUE) < 100) return(NA)  
        cor(above_bin[[gene]][nonzero_cells], above_bin$het[nonzero_cells], 
            method = "spearman", use = "pairwise.complete.obs")
      })
      
      cor_below <- sapply(gene_names, function(gene) {
        nonzero_cells <- below_bin[[gene]] > 0  
        if (sum(nonzero_cells, na.rm = TRUE) < 100) return(NA)
        cor(below_bin[[gene]][nonzero_cells], below_bin$het[nonzero_cells], 
            method = "spearman", use = "pairwise.complete.obs")
      })
      
      gene_results <- data.frame(
        bin = bin_value,
        gene = gene_names,
        corr_upper = cor_above,
        corr_below = cor_below,
        delta = cor_above - cor_below
      )
      
      results <- rbind(results, gene_results)
     
    }
  }
  results <- subset(results, !is.na(delta))
  results$corrup_arc <- (atanh(results$corr_upper))
  results$corrdown_arc <- (atanh(results$corr_below))
  
  results$delta_new <- abs(results$corrdown_arc - results$corrup_arc)
  return(results)
}

be1_delta <- correlate_het_by_bins(select(BE1, -c("cluster", "rows", "cov", "lib")))
be2_delta <- correlate_het_by_bins(select(BE2, -c("cluster", "rows", "cov", "lib")))

p1 <- ggplot(be1_delta, aes(x = as.factor(bin), y = (delta_new))) + 
  geom_boxplot(outlier.shape = NA) + 
  rasterise(geom_quasirandom(alpha = 0.01, size = 1, stroke = NA), dpi = 300) + 
  theme_classic(base_size =7, base_family = "ArialMT") + 
  ylab(expression(Delta~`correlation coefficient`)) +
  labs(x = "Heteroplasmy bin") + 
  ggtitle("143b Rep 1") +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold"))

p2 <- ggplot(be2_delta, aes(x = as.factor(bin), y = (delta_new))) + 
  geom_boxplot(outlier.shape = NA) + 
  rasterise(geom_quasirandom(alpha = 0.01, size = 1, stroke = NA), dpi = 300) + 
  theme_classic(base_size =7, base_family = "ArialMT") + 
  ylab(expression(Delta~`correlation coefficient`)) +
  labs(x = "Heteroplasmy bin") + 
  ggtitle("143b Rep 2") +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold"))

ggsave(p1, file = "~/Desktop/reznik/rrna-hotspots/revisions/results//delta_corrs_het_143b_rep1.pdf", width = 2.75, height = 2, useDingbats = FALSE, dpi = 200)
ggsave(p2, file = "~/Desktop/reznik/rrna-hotspots/revisions/results//delta_corrs_het_143b_rep2.pdf", width = 2.75, height = 2, useDingbats = FALSE, dpi = 200)