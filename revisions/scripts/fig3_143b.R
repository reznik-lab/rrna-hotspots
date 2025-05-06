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

BE1_ranks                    <- correlate_het(select(BE1, -c("cluster", "rows", "cov", "lib")))
BE2_ranks                    <- correlate_het(select(BE2, -c("cluster", "rows", "cov", "lib")))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run correlations and plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


BE1_ranksdf        <- as.data.frame(BE1_ranks)
BE1_ranksdf$gene   <- rownames(BE1_ranksdf)
BE1_ranksdf$z      <- atanh(BE1_ranksdf$BE1_ranks)

BE2_ranksdf        <- as.data.frame(BE2_ranks)
BE2_ranksdf$gene   <- rownames(BE2_ranksdf)
BE2_ranksdf$z      <- atanh(BE2_ranksdf$BE2_ranks)

BE1_ranksdf <- BE1_ranksdf[order(BE1_ranksdf$z, decreasing = TRUE),]
BE1_ranksdf$genetype <- ifelse(grepl("MT-", BE1_ranksdf$gene), "Mitochondrial", "Nuclear")
BE1_ranksdf$ranks_sort <- rank(BE1_ranksdf$BE1_ranks, ties.method = "random")

BE2_ranksdf <- BE2_ranksdf[order(BE2_ranksdf$z, decreasing = TRUE),]
BE2_ranksdf$genetype <- ifelse(grepl("MT-", BE2_ranksdf$gene), "Mitochondrial", "Nuclear")
BE2_ranksdf$ranks_sort <- rank(BE2_ranksdf$BE2_ranks, ties.method = "random")

write.csv(BE1_ranksdf, file = "~/Desktop/reznik/rrna-hotspots/revisions/results/be1_ranks.csv", row.names = FALSE)
write.csv(BE2_ranksdf, file = "~/Desktop/reznik/rrna-hotspots/revisions/results/be2_ranks.csv", row.names = FALSE)

BE1_fig3e                <- ggplot(BE1_ranksdf, aes(x = ranks_sort, y = z, colour = genetype)) + 
  geom_point(data = subset(BE1_ranksdf, genetype == "Nuclear"), size = 0.25) +
  geom_point(data = subset(BE1_ranksdf, genetype == "Mitochondrial"), alpha = 1, size = 0.5) + 
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
  ggtitle("143b Rep 1") +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = 'dashed', colour = "grey50") + 
  guides(color = guide_legend(override.aes = list(size = 1))) +
  geom_text_repel(aes(label = ifelse(gene %in% c("MT-CO1", "MT-ND2", "MT-ATP6", "MT-CYB", "RASL10A", "DCDC2B"), gene, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) +
  scale_color_manual(values = c("#E15759", "black"))

BE2_fig3e                <- ggplot(BE2_ranksdf, aes(x = ranks_sort, y = z, colour = genetype)) + 
  geom_point(data = subset(BE2_ranksdf, genetype == "Nuclear"), size = 0.25) +
  geom_point(data = subset(BE2_ranksdf, genetype == "Mitochondrial"), alpha = 1, size = 0.5) + 
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
  ggtitle("143b Rep 2") +
  geom_hline(yintercept = 0, linewidth = 0.1, linetype = 'dashed', colour = "grey50") + 
  guides(color = guide_legend(override.aes = list(size = 1))) +
  geom_text_repel(aes(label = ifelse(gene %in% c("MT-CO1", "MT-CO3", "MT-ATP6", "MT-CYB", "EBLN2", 'MYO1F'), gene, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) +
  scale_color_manual(values = c("#E15759", "black"))
ggsave(BE1_fig3e, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/be1_spearmanranks.pdf", width = 2, height = 2.75, useDingbats = FALSE, dpi = 300)

ggsave(BE2_fig3e, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/be2_spearmanranks.pdf", width = 2, height = 2.5, useDingbats = FALSE, dpi = 300)


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

BE1ranks_for_gsea <- BE1_ranksdf$z
BE2ranks_for_gsea <- BE2_ranksdf$z

names(BE1ranks_for_gsea) <- BE1_ranksdf$gene
names(BE2ranks_for_gsea) <- BE2_ranksdf$gene

BE1_ranks_for_gsea <- sort(BE1ranks_for_gsea, decreasing = TRUE)
BE2_ranks_for_gsea <- sort(BE2ranks_for_gsea, decreasing = TRUE)

BE1_gsea <- gsea(gs, BE1_ranks_for_gsea)
BE1_gsea_df <- as.data.frame(BE1_gsea)
BE1_gsea_df$padj <- as.numeric(BE1_gsea_df$padj)
BE1_gsea_df$NES <- as.numeric(BE1_gsea_df$NES)
#BE1_gsea_df <- subset(BE1_gsea_df, padj < 0.05)
BE1_gsea_df$pathway <- gsub("_", " ", BE1_gsea_df$pathway)
BE1_gsea_df$pathway <- str_to_title(BE1_gsea_df$pathway)
BE1_gsea_df$pathway <- str_sub(BE1_gsea_df$pathway, 5, nchar(BE1_gsea_df$pathway))

BE2_gsea <- gsea(gs, BE2_ranks_for_gsea)
BE2_gsea_df <- as.data.frame(BE2_gsea)
BE2_gsea_df$padj <- as.numeric(BE2_gsea_df$padj)
BE2_gsea_df$NES <- as.numeric(BE2_gsea_df$NES)
#BE2_gsea_df <- subset(BE2_gsea_df, padj < 0.05)
BE2_gsea_df$pathway <- gsub("_", " ", BE2_gsea_df$pathway)
BE2_gsea_df$pathway <- str_to_title(BE2_gsea_df$pathway)
BE2_gsea_df$pathway <- str_sub(BE2_gsea_df$pathway, 5, nchar(BE2_gsea_df$pathway))

BE1_gsea_df_subtopsignpositive <- BE1_gsea_df %>% subset(NES > 0.5) %>% arrange(padj) %>% slice_head(n =2) %>% arrange(desc(NES)) %>% mutate(color = "neg")
BE1_gsea_df_subtopsignnegative <- BE1_gsea_df %>% subset(NES < 0.5) %>% arrange(padj) %>% slice_head(n = 13) %>% arrange(desc(NES)) %>% mutate(color = "pos")
BE1_gsea_df_subtopsignpositive$pathway <- factor(BE1_gsea_df_subtopsignpositive$pathway, levels = BE1_gsea_df_subtopsignpositive$pathway)
BE1_gsea_df_subtopsignnegative$pathway <- factor(BE1_gsea_df_subtopsignnegative$pathway, levels = BE1_gsea_df_subtopsignnegative$pathway)
BE1_combined <- rbind(BE1_gsea_df_subtopsignpositive, BE1_gsea_df_subtopsignnegative)


BE2_gsea_df_subtopsignpositive <- BE2_gsea_df %>% subset(NES > 0.5) %>% arrange(padj) %>% slice_head(n =2) %>% arrange(desc(NES)) %>% mutate(color = "neg")
BE2_gsea_df_subtopsignnegative <- BE2_gsea_df %>% subset(NES < 0.5) %>% arrange(padj) %>% slice_head(n = 13) %>% arrange(desc(NES)) %>% mutate(color = "pos")
BE2_gsea_df_subtopsignpositive$pathway <- factor(BE2_gsea_df_subtopsignpositive$pathway, levels = BE2_gsea_df_subtopsignpositive$pathway)
BE2_gsea_df_subtopsignnegative$pathway <- factor(BE2_gsea_df_subtopsignnegative$pathway, levels = BE2_gsea_df_subtopsignnegative$pathway)
BE2_combined <- rbind(BE2_gsea_df_subtopsignpositive, BE2_gsea_df_subtopsignnegative)

BE1_gsea_plot <- ggplot(BE1_combined, aes(x = pathway, y = NES, fill = color)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme_classic(base_size = 7, base_family = "ArialMT") + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        axis.ticks.y = element_blank(),  
        legend.position = "none",
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("143b Rep 1") +
  labs(x = "", y = "Normalized enrichment score") + 
  scale_fill_manual(values = c("#4E79A7", "#E15759"))


BE2_gsea_plot <- ggplot(BE2_combined, aes(x = pathway, y = NES, fill = color)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme_classic(base_size = 7, base_family = "ArialMT") + 
  theme(axis.line = element_line(linewidth = 0.1),
        axis.ticks = element_line(linewidth = 0.05),
        axis.ticks.length = unit(1, "pt"),
        axis.ticks.y = element_blank(),  
        legend.position = "none",
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold")) + 
  ggtitle("143b Rep 2") +
  labs(x = "", y = "Normalized enrichment score") + 
  scale_fill_manual(values = c("#4E79A7", "#E15759"))

ggsave(BE1_gsea_plot, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/BE1_gsea.pdf", width =3.5, height = 2.5)
ggsave(BE2_gsea_plot, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/BE2_gsea.pdf", width =3.5, height = 2.5)

BE1_gsea_df$colour <- ifelse(BE1_gsea_df$padj < 0.05 & BE1_gsea_df$NES > 1, "pos", 
                             ifelse(BE1_gsea_df$padj < 0.05 & BE1_gsea_df$NES < -1, "neg", "nothing"))
be1_gsea <- ggplot(BE1_gsea_df, aes(y = -log10(padj), x = NES, colour = colour)) + 
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
        legend.position = "none") +
  scale_colour_manual(values = c("#4E79A7", "black", "#E15759" )) + 
  geom_text_repel(aes(label = ifelse(pathway %in% c(" Mitochondrial Small Ribosomal Subunit", 
                                                    " Respiratory Chain Complex Iv",
                                                    " Nmda Glutamate Receptor Activity",
                                                    " Gated Channel Activity",
                                                    " Proton Transporting Atp Synthase Complex"), pathway, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) 

ggsave(be1_gsea, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/be1_gsea_volcano.pdf", width = 3.5, height = 2.5, useDingbats = FALSE, dpi = 300)

BE2_gsea_df$colour <- ifelse(BE2_gsea_df$padj < 0.05 & BE2_gsea_df$NES > 1, "pos", 
                             ifelse(BE2_gsea_df$padj < 0.05 & BE2_gsea_df$NES < -1, "neg", "nothing"))
be2_gsea <- ggplot(BE2_gsea_df, aes(y = -log10(padj), x = NES, colour = colour)) + 
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
        plot.title = element_text(size = 6, family = "ArialMT", hjust = 0.5, face = "bold"),
        legend.position = "none") +
  scale_colour_manual(values = c("#4E79A7", "black", "#E15759" )) + 
  ggtitle("143b Rep 2") +
  geom_text_repel(aes(label = ifelse(pathway %in% c(" Mitochondrial Small Ribosomal Subunit", 
                                                    " Respiratory Chain Complex Iv",
                                                    " Mitochondrial Electron Transport Cytochrome C To Oxygen",
                                                    " Proton Transporting Atp Synthase Complex"), pathway, "")), family = 'ArialMT', size = 2, max.overlaps = 10000, segment.size = 0.1) 

ggsave(be2_gsea, file = "~/Desktop/reznik/rrna-hotspots/revisions//figure_panels/fig3/be2_gsea_volcano.pdf", width = 3, height = 2, useDingbats = FALSE, dpi = 300)
