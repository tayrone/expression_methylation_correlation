library(tidyverse)
library(gdata)


load("./rdatas_methylation/g3_signature.RData")


hits_expression <- pheno_ids %>% inner_join(tibble(gene = hits)) %>% 
  select(expression = exp_values, everything())


gdata::keep(hits_expression, sure = T)


load("./rdatas_methylation/group3_probewised.RData")

methylation_coefficient <- as.data.frame(fit[["coefficients"]]) %>% 
  rownames_to_column(var = "probe_id") %>% 
  select(everything(), methylation = `control-Group3`)
  

test <- diff_methylated %>% 
  as_tibble() %>% 
  select(probe_id = name, probe_location = ucsc_refgene_group, 
         gene = ucsc_refgene_name, aveexpr) %>% 
  inner_join(methylation_coefficient)

hist(test$methylation)
hist(test$aveexpr)
