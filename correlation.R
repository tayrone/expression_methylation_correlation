library(tidyverse)

#---- Loads and organizes expression data for further use ----

load("./rdatas/shh_network_dictionary.RData")

colnames(gexp) <- str_remove(colnames(gexp), ".CEL")

dictionary <- dictionary %>% mutate(across(probe_id, as.character))

gexp <- gexp %>% 
  rownames_to_column(var = "probe_id") %>% 
  inner_join(dictionary) %>% 
  select(probe_id, gene = gene_symbol, everything())

gdata::keep(gexp, sure = T)


#---- Organizes methylation data for further use ----

load("./rdatas/shh_filtered.RData")
load("./rdatas_methylation/shh_probewised.RData")

beta_values <- as.data.frame(getBeta(filtered_quantile_normalized))

diff_methylated <- diff_methylated %>% 
  as_tibble() %>% 
  select(probe_id = name, probe_location = ucsc_refgene_group, 
         gene = ucsc_refgene_name)

colnames(beta_values) <- str_extract(colnames(beta_values), "[^_]+")

beta_values <- beta_values %>% 
  rownames_to_column(var = "probe_id") %>%
  inner_join(diff_methylated) %>% 
  select(probe_id, probe_location, gene, everything()) %>% 
  drop_na(gene)

gdata::keep(gexp, beta_values, sure = T)


#---- 
# Implementar teste de correlação em si. Talvez seja interessante aumentar o
# data frame gexp para que ele possa ser pareado com o beta_values, linha a linha.

cor.test(beta_values[, -c("probe_id", "probe_location", "gene")])

corr = cor.test(as.numeric(bVals_Tumor[which(rownames(bVals_Tumor) == bValsPlot$sonda[i]),]), 
                as.numeric(signature_gexp[which(signature_dictionary$gene_symbol == bValsPlot$gene[i]),]),
                method = "spearman")