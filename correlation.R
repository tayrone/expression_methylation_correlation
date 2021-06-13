library(tidyverse)
library(gdata)
library(GEOquery)

# Set subgroup first
investigated_subgroup <- c("Group4" = "g4")

load(paste0("../rdatas/", investigated_subgroup, "_network_dictionary.RData"))
load(paste0("../rdatas/", investigated_subgroup, "_filtered.RData"))
load(paste0("../rdatas/", investigated_subgroup, "_probewised.RData"))


#---- Loads and organizes expression data for further use ----

colnames(gexp) <- str_remove(colnames(gexp), ".CEL")

dictionary <- dictionary %>% mutate(across(probe_id, as.character))

gexp <- gexp %>% 
  rownames_to_column(var = "probe_id") %>% 
  inner_join(dictionary) %>% 
  select(exp_probe_id = probe_id, gene = gene_symbol, everything())


#---- Organizes methylation data for further use ----

beta_values <- as.data.frame(getBeta(filtered_quantile_normalized))

diff_methylated <- diff_methylated %>% 
  as_tibble() %>% 
  select(probe_id = name, probe_location = ucsc_refgene_group, 
         gene = ucsc_refgene_name)

colnames(beta_values) <- str_extract(colnames(beta_values), "[^_]+")

beta_values <- beta_values %>% 
  rownames_to_column(var = "probe_id") %>%
  inner_join(diff_methylated) %>% 
  select(met_probe_id = probe_id, probe_location, gene, everything()) %>% 
  drop_na(gene)


#---- Prepares expression and beta values data frame for correlation test ----

genes_gexp <- select(gexp, exp_probe_id, gene)

met_genes_without_exp_probe <- beta_values %>% 
  select(met_probe_id, probe_location, gene) %>% 
  left_join(genes_gexp) %>% 
  filter(is.na(exp_probe_id)) %>% 
  select(gene) %>% 
  distinct() %>% 
  pull(gene)

beta_values <- beta_values %>% 
  filter(!(gene %in% met_genes_without_exp_probe)) %>% 
  arrange(gene)


genes_met <- select(beta_values, met_probe_id, gene)

exp_genes_without_met_probe <- gexp %>% 
  select(exp_probe_id, gene) %>% 
  left_join(genes_met) %>% 
  filter(is.na(met_probe_id)) %>% 
  select(gene) %>% 
  distinct() %>% 
  pull(gene)

gexp <- gexp %>% 
  filter(!(gene %in% exp_genes_without_met_probe))


beta_gene_format <- beta_values %>% 
  select(gene) 

gexp <- inner_join(gexp, beta_gene_format) %>% 
  arrange(gene)

gdata::keep(gexp, beta_values, investigated_subgroup, sure = T)


#---- Identifies which samples did not pass the quality control on either
# methylation and expressions analyses, so all of them can be discarded ----

accession_expression <- getGEO("GSE85217", GSEMatrix = TRUE, getGPL = FALSE)

accession_expression <- 
  accession_expression[["GSE85217_series_matrix.txt.gz"]]@phenoData@data %>% 
  select(title, exp_accession = geo_accession, subgroup = `subgroup:ch1`,
         subtype = `subtype:ch1`) %>% 
  filter(subgroup == names(investigated_subgroup)) %>% 
  as_tibble()


accession_methylation <- read_csv("./csvs/accession_methylation.csv") %>% 
  select(met_accession = Accession, title = Title)


accession_dictionary <- accession_methylation %>% 
  mutate(title = str_remove(title, "_methylation")) %>% 
  inner_join(accession_expression) %>% 
  select(title, everything())

gdata::keep(accession_dictionary, beta_values, gexp, sure = T)


beta_gsms <- beta_values %>% 
  select(starts_with("GSM")) %>% 
  colnames

gexp_gsms <- gexp %>% 
  select(starts_with("GSM")) %>% 
  colnames


# Removes control samples attached to case samples on previous steps.
# Gexp has no control samples, since controls are taken into account in
# the signature definition.

met_controls <- which(!(beta_gsms %in% accession_dictionary$met_accession))
beta_gsms <- beta_gsms[-met_controls]
rm(met_controls)


accession_dictionary <- accession_dictionary %>% 
  filter(exp_accession %in% gexp_gsms) %>% 
  filter(met_accession %in% beta_gsms)


#---- Keeps only samples present in the dictionary, while also
# removing annotating columns ---- 

beta_annotation <- beta_values %>% 
  select(!starts_with("GSM"))

gexp_annotation <- gexp %>% 
  select(!starts_with("GSM"))


beta_values <- beta_values %>% 
  select(all_of(pull(accession_dictionary, met_accession)))

gexp <- gexp %>% 
  select(all_of(pull(accession_dictionary, exp_accession))) 


gdata::keep(accession_dictionary, beta_values, gexp,  
            beta_annotation, gexp_annotation, sure = T)

#---- Correlation values are finally calculated, row by row, taking both 
# tibbles in parallel ----

indexes <- 1:nrow(gexp)

corr <- function(row){
  cor.test(as.numeric(dplyr::slice(beta_values, row)), 
           as.numeric(dplyr::slice(gexp, row)), 
           method = "spearman", exact = F) 
}

cor_result <- map(indexes, corr) %>% 
  map(.f = function(value) value[["estimate"]]) %>% 
  unlist()

beta_annotation <- beta_annotation %>% 
  mutate(correlation = cor_result)

