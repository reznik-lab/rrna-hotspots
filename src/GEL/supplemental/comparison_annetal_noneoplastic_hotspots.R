## AUTHOR: SONIA BOSCENCO

indel_hotspots        <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/supp_tables/table_3_indel_hotspots.txt", sep = " ")
snv_hotspots          <- read.csv("~/Desktop/reznik/rrna-hotspots/revisions/data/hotspots_snv_5.txt", sep = "\t")

juetalhots            <- readxl::read_excel("~/Downloads/41588_2024_1838_MOESM4_ESM.xlsx", sheet = 6,skip = 3)
juetalhots            <- juetalhots %>% subset(`Dloop                    (0=no, 1=yes)` == 0)

juetalhots$mutrate    <- indel_hotspots$mutant[match(juetalhots$POS, indel_hotspots$start)]
juetalhots$mutrate    <- ifelse(is.na(juetalhots$mutrate), snv_hotspots$mutations_at_pos[match(juetalhots$POS, snv_hotspots$pos)], juetalhots$mutrate)

ggplot(juetalhots, aes(x = mutrate, y = patientN)) + 
  geom_point()

juetalhots$isours <- ifelse(!is.na(juetalhots$mutrate), TRUE, FALSE)

tbl_rate <- as.data.frame(table(juetalhots$isours))

snv_hotspots$inju <- ifelse(snv_hotspots$pos %in% juetalhots$POS, TRUE, FALSE)
sum(snv_hotspots$inju)

indel_hotspots$inju <- ifelse(indel_hotspots$start %in% juetalhots$POS, TRUE, FALSE)
sum(indel_hotspots$inju)