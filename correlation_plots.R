library(tidyverse)
library(GEOquery)
library(gdata)

investigated_subgroup <- c("Group4" = "g4")
plots_pallete <- c("#008AB8", "#C70039")


interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

load(paste0("../rdatas/", investigated_subgroup, "_probewised.RData"))

load(paste0("./correlation_results/", investigated_subgroup, "_results.RData"))

gdata::keep(beta_annotation, plots_pallete, investigated_subgroup, interesse_regs,
            sure = T)


#---- Minor alterations on tibble to be manipulated ----

theme_set(theme_light())

beta_annotation <- beta_annotation %>% 
  mutate(probe_location = str_extract(probe_location, "[^;]+"),
         regulation = ifelse(correlation > 0, "up", "down"))


#---- Number of probes by location ----

first_plot <- beta_annotation %>% 
  group_by(probe_location, regulation) %>% 
  summarise(n_probes = n())
   
first_plot %>% 
  ggplot() +
  geom_col(aes(probe_location, n_probes, fill = regulation), alpha  = 0.85, 
           position = "dodge") +
  labs(x = "", y = "", fill = "Regulation", 
       title = "Amount of methylation probes by location") +
  scale_fill_manual(values = plots_pallete)

  ggsave(paste0("./plots/", investigated_subgroup, "/probes_count.png"),
         width = 12, height = 8)


#---- Overall count of up and down regulated probes ----

beta_annotation %>% 
  group_by(gene) %>% 
  summarise(regulation_summary = ifelse(n_distinct(regulation) == 2, "both", "one"),
            regulation_summary = ifelse(regulation_summary == "one" & regulation == "up", 
                                        "up", regulation_summary),
            regulation_summary = ifelse(regulation_summary == "one" & regulation == "down",
                                        "down", regulation_summary)) %>% 
  ggplot() +
    geom_histogram(aes(regulation_summary), stat = "count") 


#---- Number of probes by location, for regulons of interest ----

second_plot <- beta_annotation %>% 
  filter(gene %in% interesse_regs) %>% 
  group_by(probe_location, regulation) %>% 
  summarise(n_probes = n())

second_plot %>% 
  ggplot() +
  geom_col(aes(probe_location, n_probes, fill = regulation), alpha  = 0.85, 
           position = "dodge") +
  labs(x = "", y = "", fill = "Regulation", 
       title = "Methylation probes by location, for regulons of interest") +
  scale_fill_manual(values = plots_pallete)

ggsave(paste0("./plots/", investigated_subgroup, "/probe_count_regs_interest.png"),
       width = 12, height = 8)
  

#---- Distribution of correlation values for each regulon of interest ---

beta_annotation %>% 
  filter(gene %in% interesse_regs) %>% 
  group_by(gene) %>% 
  ggplot(aes(gene, correlation)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "", y = "Correlation Coefficient",
       title = "Distribution of correlation values for each regulon of interest")

ggsave(paste0("./plots/", investigated_subgroup, "/correlation_distribution.png"),
       width = 12, height = 8)




