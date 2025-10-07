### AUGUST 13th 2024
## AUTHOR: SONIA BOSCENCO
## SCRIPT TO ANALYZE RNA-SEQ DATA FROM HEK293 CELLS
so <- readRDS("~/Desktop/reznik/rrna-hotspots/hek_cells/data/hek_processed_seurat.rds")
library(ggpubr)
source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')
library(fgsea)
library(qusage)
library(Seurat)
library(ggrastr)
library(ggbeeswarm)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot and save UMAP FIG 3B
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig3b                <- DimPlot(so, group.by = "seurat_clusters", label = TRUE, label.size = 3) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  theme_classic(base_size = 7, 
                base_family = 'ArialMT') +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed", length = unit(8,'pt'))))

if(save) ggsave(fig3b, file = "~/Desktop/reznik/1227_atac/figures/fig3b_umap.pdf", width = 2.5, height = 2.5)
rm(fig3b, gex_mat, het_df)

# looking at heteroplasmy distribution, clustering is not based on this 
FeaturePlot(so, "af1227") + scale_color_viridis() + theme_classic(base_family = "ArialMT", base_size = 7) + labs(title = "") + 
  theme(axis.line = element_line(linewidth = 0.1),
        legend.direction = "horizontal",
        legend.ticks = element_blank(),
        legend.key.height = unit(6, "pt"),
        legend.key.width = unit(12, "pt"),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        legend.position = "bottom") 

all.genes            <- rownames(gex_mat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to correlate heteroplasmy to gene expression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

colnames(het)            <- c("het")
colnames(cluster)        <- c("cluster")
colnames(cov)            <- c("cov")

gene_expression_mat$rows <- rownames(gene_expression_mat)
het$rows                 <- rownames(het)
cluster$rows             <- rownames(cluster)
cov$rows                 <- rownames(cov)

master_df                <- merge(gene_expression_mat, het, by = "rows")
master_df                <- merge(master_df, cluster, by = "rows")
master_df                <- merge(master_df, cov, by = "rows")

het_master_1227          <- master_df %>% subset(het > 0.05)
ranks                    <- correlate_het(select(het_master_1227, -c("cluster", "rows", "cov")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run correlations and plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


ranksdf        <- as.data.frame(ranks)
ranksdf$gene   <- rownames(ranksdf)
ranksdf$z      <- atanh(ranksdf$ranks)


write.csv(ranksdf, file = "~/Desktop/reznik/rrna-hotspots/revisions/results/hek_ranks.csv", row.names = FALSE)


p <- ggplot(subset(het_master_1227, `MT-CO1` > 0), aes(x = het, y = `MT-CO1`)) + 
  geom_point(alpha = 0.5, colour = "#B07AA1", stroke = NA, size = 1) + 
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
  ggtitle("HEK293") 


ggsave(p, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/fig3b_heteroplasmy_hek293.pdf", width = 2, height = 1.25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot ranked spearman correlations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
ranksdf <- ranksdf[order(ranksdf$z, decreasing = TRUE),]
ranksdf$genetype <- ifelse(grepl("MT-", ranksdf$gene), "Mitochondrial", "Nuclear")

ranksdf$ranks_sort <- rank(ranksdf$ranks, ties.method = "random")

fig3e                <- ggplot(ranksdf, aes(x = ranks_sort, y = z, colour = genetype)) + 
  geom_point(data = subset(ranksdf, genetype == "Nuclear"), size = 0.25) +
  geom_point(data = subset(ranksdf, genetype == "Mitochondrial"), alpha = 1, size = 0.5) + 
  labs(x = "", y = "Fisher transformed spearman coefficient", colour = "")  + 
  theme_classic(base_size = 7, 
                base_family = 'ArialMT') + 
  theme(legend.position = "bottom",
        legend.key.size = unit(6, "pt"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("HEK293") +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = 'dashed', colour = "grey50") + 
  guides(color = guide_legend(override.aes = list(size = 1))) +
  geom_text_repel(aes(label = ifelse(gene %in% c("MT-CO1", "MT-ND5", "MT-ND2", "MT-ATP6", "MT-CYB", "CAMK2N2", "DCDC2B"), gene, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) +
  scale_color_manual(values = c("#E15759", "black"))

ggsave(fig3e, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/hek_293_spearmanranks.pdf", width = 2, height = 2.5, useDingbats = FALSE, dpi = 300)

gmtfile                  <- "~/Desktop/reznik/reference/gsea_terms/c5.go.v2023.1.Hs.symbols.gmt"
gs                       <-  read.gmt(gmtfile)
gsea                   <- function(gs, ranks){
  resgsea = fgsea(pathways = gs,stats = ranks,minSize=1,maxSize=500)
  resgsea = resgsea[order(resgsea$padj,decreasing = FALSE),]
  resdf = apply(resgsea,2,as.character)
  resdf[,'padj'] <- as.numeric(resdf[,'padj'])
  resdf[,'NES'] <- as.numeric(resdf[,'NES'])
  #resdf_sig = resdf[which(resdf[,'padj'] < .05),]
  return(resdf)
}

ranks_for_gsea <- ranksdf$z

names(ranks_for_gsea) <- ranksdf$gene
ranks_for_gsea <- sort(ranks_for_gsea, decreasing = TRUE)
gsea <- gsea(gs, ranks_for_gsea)
gsea_df <- as.data.frame(gsea)
gsea_df$padj <- as.numeric(gsea_df$padj)
gsea_df$NES <- as.numeric(gsea_df$NES)
#gsea_df <- subset(gsea_df, padj < 0.05)
gsea_df$pathway <- gsub("_", " ", gsea_df$pathway)
gsea_df$pathway <- str_to_title(gsea_df$pathway)
gsea_df$pathway <- str_sub(gsea_df$pathway, 5, nchar(gsea_df$pathway))


gsea_df$colour <- ifelse(gsea_df$padj < 0.05 & gsea_df$NES > 1, "pos", 
                             ifelse(gsea_df$padj < 0.05 & gsea_df$NES < -1, "neg", "nothing"))
hek_gsea <- ggplot(gsea_df, aes(y = -log10(padj), x = NES, colour = colour)) + 
  geom_point(size = 1, alpha = 0.5, stroke = NA) + 
  geom_hline(yintercept = -log10(0.05), linewidth = 0.1, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1, linewidth = 0.1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -1, linewidth = 0.1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.1, linetype = "dashed", colour = "grey50") +
  theme_classic(base_size = 7, base_family = "ArialMT") + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        legend.position = "none",
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) +
  scale_colour_manual(values = c("#4E79A7", "black", "#E15759" )) + 
  ggtitle("HEK293") +
  geom_text_repel(aes(label = ifelse(pathway %in% c(" Cytosolic Small Ribosomal Subunit", 
                                                    " Respiratory Chain Complex Iv",
                                                    " Axonemal Dynein Complex",
                                                    " Proton Transporting Atp Synthase Complex"), pathway, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) 


ggsave(hek_gsea, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/hek_gsea.pdf", width =3, height = 2)

write.csv(gsea_df, file = "~/Desktop/reznik/rrna-hotspots/hek_cells/results/go_gsea.csv", row.names = FALSE)

gsea_df_subtopsignpositive <- gsea_df %>% subset(NES > 0.5) %>% arrange(padj) %>% slice_head(n =2) %>% arrange(desc(NES)) %>% mutate(color = "neg")

gsea_df_subtopsignnegative <- gsea_df %>% subset(NES < 0.5) %>% arrange(padj) %>% slice_head(n = 13) %>% arrange(desc(NES)) %>% mutate(color = "pos")

gsea_df_subtopsignpositive$pathway <- factor(gsea_df_subtopsignpositive$pathway, levels = gsea_df_subtopsignpositive$pathway)

gsea_df_subtopsignnegative$pathway <- factor(gsea_df_subtopsignnegative$pathway, levels = gsea_df_subtopsignnegative$pathway)

combined <- rbind(gsea_df_subtopsignpositive, gsea_df_subtopsignnegative)

gsea_plot <- ggplot(combined, aes(x = pathway, y = NES, fill = color)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
 theme_classic(base_size = 7, base_family = "ArialMT") + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        axis.ticks.y = element_blank(),  
        legend.position = "none",
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("HEK293") +
   labs(x = "", y = "Normalized enrichment score") + 
  scale_fill_manual(values = c("#4E79A7", "#E15759"))



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
  
  return(results)
}
result_df <- correlate_het_by_bins(select(het_master_1227, -c("cluster", "rows", "cov")))
result_df <- subset(result_df, !is.na(delta))
result_df$corrup_arc <- (atanh(result_df$corr_upper))
result_df$corrdown_arc <- (atanh(result_df$corr_below))

result_df$delta_new <- abs(result_df$corrdown_arc - result_df$corrup_arc)
p1 <- ggplot(result_df, aes(x = as.factor(bin), y = (delta_new))) + 
  geom_boxplot(outlier.shape = NA) + 
  rasterise(geom_quasirandom(alpha = 0.01, size = 1, stroke = NA), dpi = 300) + 
  theme_classic(base_family = "ArialMT", base_size = 7) + 
  ylab(expression(Delta~`correlation coefficient`)) +
labs(x = "Heteroplasmy bin") + 
  ggtitle("HEK293") +
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold"))

ggsave(p1, file = "~/Desktop/reznik/rrna-hotspots/revisions/results/hek293_heteroplasmythresholdbins.pdf", width = 2.75, height = 2, useDingbats = FALSE, dpi = 200)


