rm(list = ls(all.names = TRUE))
gc()

.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)
library(fontawesome)
library(data.table)

## input files
mutation_file        <- '~/re_gecip/cancer_pan/boscens/data/filtervars_0223.tsv'
fa_file              <- "~/re_gecip/cancer_pan/boscens/data/GRCh38_chrM.fa"
chrM_annotated_file  <- '~/re_gecip/cancer_pan/boscens/data/chrM_annotated.txt'

## output files
hotspots_output      <- '~/re_gecip/cancer_pan/boscens/data/hotspots_snv_0220.txt' 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tabulate the 64 unique possible trinucleotides for all retained positions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fa                   <- paste(readLines(fa_file), collapse = '')
dat                  <- data.table(pos=1:16569)
dat$nt               <- strsplit(fa,'')[[1]]
dat$ntMinus1         <- c('-',dat$nt[1:(nrow(dat)-1)])
dat$ntPlus1          <- c(dat$nt[2:(nrow(dat))],'-')
dat$trinuc           <- paste0(dat$ntMinus1,dat$nt,dat$ntPlus1)
dat                  <- dat[,c('pos','trinuc'),with=F]

d                    <- fread(chrM_annotated_file)
d                    <- merge(d, dat, by='pos', all.x=T)
d$trinuc             <- as.character(d$trinuc)
d$pos                <- as.integer(d$pos )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate trinucleotide mutabilities
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maf                  <- read.csv(mutation_file, sep = "\t", header = T)
maf$Start_Position   <- as.integer(maf$Start_Position)
maf$flanking_bps     <- d$trinuc[match(maf$Start_Position, d$pos)]
maf$Hugo_Symbol      <- d$symbol[match(maf$Start_Position, d$pos)]
maf                  <- maf %>% 
                        subset(Variant_Type == "SNP")

tbl                  <- as.data.frame(table(maf$flanking_bps))
colnames(tbl)        <- c('trinuc','s_trinuc')
s_total              <- sum(tbl$s_trinuc)
tbl$trinuc           <- as.character(tbl$trinuc)
d$s_trinuc           <- tbl$s_trinuc[match(d$trinuc, tbl$trinuc)]
d[is.na(s_trinuc),]  <- 0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add number of samples mutated per gene
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get mutations per position
mutation_counts      <- as.data.frame(table(maf$Start_Position))
d$mutations_at_pos   <- mutation_counts$Freq[match(d$pos, mutation_counts$Var1)]
d$mutations_at_pos[is.na(d$mutations_at_pos)]<- 0



## add mutations per gene
tbl                  <- as.data.frame(table(maf$Hugo_Symbol))
setnames(tbl,'Freq','mutations_in_gene')
setnames(tbl, "Var1", "symbol")
d$mutations_in_gene  <- tbl$mutations_in_gene[match(d$symbol, tbl$symbol)]
d$mutations_in_gene[is.na(d$mutations_in_gene)] <- 0
d$mu_pos             <- d$s_trinuc / s_total


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# normalize the position mutability by the gene's total trinucleotide-mutability
# this is the binomial parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- d %>% subset(mutations_at_pos > 0)
get_binomial_parameter <- function(res) {
  res$binomial_probability <- res$mu_pos / sum(res$mu_pos)
  res
}
d <- d[,get_binomial_parameter(.SD),by=symbol]
d <- d[order(d$pos),]

d <- d %>% subset(mutations_in_gene > 0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test each position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_position <- function(dd) {
  successes <- dd$mutations_at_pos
  trials <- dd$mutations_in_gene
  probability <- dd$binomial_probability
  expected_number_muts_gene <- trials*probability
  if(successes < 2*expected_number_muts_gene){
    p.value <- 1
  } else{
    p.value <- binom.test(x=successes, n=trials, p=probability, alternative="greater")$p.value
  }
  out <- list(p.value=p.value)
  out
}
hs <- d[, test_position(.SD), by=pos]
res <- merge(d, hs, by='pos', all=F)
res$q.value <- p.adjust(res$p.value,method='BH')
res <- res[order(p.value, decreasing=F),]
res_sign <- subset(res, q.value < 0.05)

write.table(res, file = "/re_gecip/cancer_pan/boscens/data/hotspots_snv_0223.txt", sep = "\t")

res_filtered <- res %>% subset(q.value < 0.05)
