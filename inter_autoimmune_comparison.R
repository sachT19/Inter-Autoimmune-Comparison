#Sachleen Tuteja, May/June 2022, inter-autoimmune comparison and controls
library(tidyverse)
library(tidyr)
library(ggsci)
library(ggpubr)
library(VennDiagram)
library(RColorBrewer)
library(VennDiagram)
library(stringr)
library(qqman)

#lupus data 
lupus = readRDS("/Users/sachleentuteja/Desktop/northwestern/prs/raw/published_prs.RDS")
list2env(lupus, envir = .GlobalEnv) 
gwas328 <- read.csv(file = '/Users/sachleentuteja/Desktop/northwestern/prs/raw/GCST90011866_buildGRCh37.tsv', 
                     sep = '\t', header = TRUE)

#ra data
ra = readRDS("/Users/sachleentuteja/Desktop/northwestern/prs/raw/ra/rheumatoid_arthritis_pgs.RDS")
list2env(ra, envir = .GlobalEnv)

#manhattan plot for gwas for PGS000328
manhattan(gwas328, chr="chromosome", bp="base_pair_location", snp = "variant_id",
          beta = "beta", p="p_value", ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, 
          colors = c("#666666", "#CC6600"), pch=20, 
          genomewideline = F, suggestiveline = F)

#sle pgs: "PGS000328" "PGS000754" "PGS000772"
#ra pgs: "PGS000194" "PGS001875" "PGS002260"

#PGS000328 - GRCh38
#PGS000194 - GRCh38
PGS000328$logweight = log(PGS000328$effect_weight)
PGS000194$logweight = log(PGS000194$effect_weight)
shared_sle_ra1 = inner_join(PGS000328, 
                           PGS000194, 
                           by  = "rsID")
myCol <- c('#000000')
ggscatter(shared_sle_ra1, x = "logweight.x", y = "logweight.y",
          add = "reg.line", xlab = "Lupus Genetic Risk Score", ylab = "RA Genetic Risk Score", 
          color = myCol)
ggscatter(shared_sle_ra1, pch = 1, x = "effect_weight.x", y = "effect_weight.y",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "Lupus Genetic Risk Score", 
          ylab = "RA Genetic Risk Score", 
          color = myCol)

#PGS000754 - GRCh37
#PGS001875 - GRCh37
PGS000754$logweight = log(PGS000754$effect_weight)
PGS001875$logweight = log(PGS001875$effect_weight)
shared_sle_ra2 = inner_join(PGS000754,
                            PGS001875,
                            by = "rsID")
ggscatter(shared_sle_ra2, x = "logweight.x", y = "logweight.y",
          add = "reg.line", xlab = "Lupus Genetic Risk Score", ylab = "RA Genetic Risk Score", 
          color = myCol)
ggscatter(shared_sle_ra2, pch = 1, x = "effect_weight.x", y = "effect_weight.y",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "Lupus Genetic Risk Score", 
          ylab = "RA Genetic Risk Score", 
          color = myCol)

#PGS000772 - GRCh37
#PGS002260 - GRCh37
PGS000772$logweight = log(PGS000772$effect_weight)
PGS002260$logweight = log(PGS002260$effect_weight)
shared_sle_ra3 = inner_join(PGS000772,
                            PGS002260,
                            by = c("chr", "chr_position"))
#joined = prs1 %>% 
  #inner_join(prs2, by = c("chr", "pos"))
ggscatter(shared_sle_ra3, x = "logweight.x", y = "logweight.y",
          add = "reg.line", xlab = "Lupus Genetic Risk Score", ylab = "RA Genetic Risk Score", 
          color = myCol)
ggscatter(shared_sle_ra3, pch = 1, x = "effect_weight.x", y = "effect_weight.y",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", 
          xlab = "Lupus Genetic Risk Score", 
          ylab = "RA Genetic Risk Score", 
          color = myCol)