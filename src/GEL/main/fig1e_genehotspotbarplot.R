## AUTHOR: SONIA BOSCENCO 
rm(list=ls())
source('~/Desktop/reznik/rrna-hotspots/revisions/scripts//prerequisites.R')

hotspots        <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t")


gene_list_1e       <- copy(gene_list)
gene_list_1e$gene  <- paste0("MT-", gene_list_1e$gene)

table_pergene      <- as.data.frame(table(hotspots$symbol))

# sort by frequency 
table_pergene      <- merge(table_pergene, gene_list_1e, by.x = "Var1", by.y = "gene")

table_pergene       <- table_pergene %>% 
                       arrange(Freq)

table_pergene$Var1  <- factor(table_pergene$Var1, levels = table_pergene$Var1)
table_pergene$color <- ifelse(table_pergene$Var1 %in% c("MT-RNR1", "MT-RNR2"), "#4E79A7", "#BAB0AC")

panel_1e            <- ggplot(table_pergene, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", aes(fill = color)) + 
  theme_std() + 
  coord_flip() + 
  scale_y_continuous(expand = c(0,0), position = "right")+ 
  scale_x_discrete(expand = c(0,0))+
  scale_fill_manual(values = c("#4E79A7", "#BAB0AC")) +
  labs(x = "", 
       y = "Number of hotspots in gene",
       fill = "") +
  theme(legend.box = element_blank(),
        legend.text = element_blank(),
        legend.key = element_blank(),
        legend.position = "none",
        axis.ticks = element_line(linewidth = 0.05),
        axis.line = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(margin = margin(r = 1)),
        axis.ticks.length = unit(1, "pt"))

ggsave(panel_1e,file = '~/Desktop/reznik/rrna-hotspots/revisions/figure_panels/panel_1e.pdf',width=2,height=2.5)