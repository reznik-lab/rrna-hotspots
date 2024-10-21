## set seed for reproducibility of figures
set.seed(123)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

required.packages <- c('data.table','ggplot2', 'dplyr', 'tidyr', 'stringr', 'RColorBrewer', 'ggsignif', 'ggrepel', 'binom', 'reshape2', 'ggpubr')
hide <- suppressMessages(lapply(required.packages, require, character.only = TRUE))
missing.packages <- required.packages[!required.packages %in% (.packages())]
if(length(missing.packages)>0) stop(paste('Could not load required packages:',paste(missing.packages,collapse=', ')))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define commonly used objects and helper functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coding_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','Frame_Shift_Ins','Translation_Start_Site','Nonstop_Mutation','Silent')
truncating_classes <- c('Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site')

## ggplot theme
theme_std <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
      line = element_line(colour = "black", linetype = 1))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define mtdna gene_list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gene_list <- read.table(text = "gene start end
1 CR 1 576
17 TF 577 647
19 RNR1 648 1601
20 TV 1602 1670
21 RNR2 1671 3229
25 TL1 3230 3304
27 ND1 3307 4262
28 TI 4263 4330
31 TM 4402 4469
32 ND2 4470 5511
33 TW 5512 5579
38 OLR 5725 5780
42 CO1 5904 7445
45 TD 7518 7585
46 CO2 7586 8269
48 TK 8295 8364
50 ATP8 8366 8550
51 ATP6 8550 9207
52 CO3 9207 9990
53 TG 9991 10058
54 ND3 10059 10404
55 TR 10405 10469
56 ND4L 10470 10760
57 ND4 10760 12137
58 TH 12138 12206
59 TS2 12207 12265
60 TL2 12266 12336
61 ND5 12337 14148
65 CYB 14747 15887
66 TT 15888 15953
2 CR 16024 16569
29 TQ 4330 4400
35 TA 5587 5655
37 TN 5657 5725
39 TC 5780 5826
40 TY 5826 5891
43 TS1 7446 7514
62 ND6 14149 14673
63 TE 14674 14742
70 TP 15956 16023", header = TRUE, sep = " ")
gene_list$length    <- gene_list$end - gene_list$start + 1

## no tRNA
filtered_gene_list  <- gene_list[!grepl("^T", gene_list$gene), ]