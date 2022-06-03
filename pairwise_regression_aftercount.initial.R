#!/usr/bin/env Rscript

if (!require("dplyr")) install.packages("dplyr")
if (!require("data.table")) install.packages("data.table")
if (!require("optparse")) install.packages("optparse")
if (!require("tidyr")) install.packages("tidyr")
if (!require("broom")) install.packages("broom")

library(optparse)
library(data.table)
library(tidyr)
library(broom)
library(dplyr)

option_list = list(
    make_option(c("-f", "--counts_file"), type="character", default=NULL, 
              help="filtered counts file name", metavar="character"),
    make_option(c("-t", "--totalsums_out"), type="character", default="pairwise_totalsums.txt",
              help="totalsums output file name [default= %default]", metavar="character")  
    make_option(c("-o", "--lm_outfile"), type="character", default="pairwise_regtests.txt", 
              help="regtest output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts_file) || is.null(opt$lm_outfile)){
    print_help(opt_parser)
    stop("At least two arguments must be supplied (input & output file).n", call.=FALSE)
}

transition_table <- data.table::fread("transition_table.txt",
                                      header=T, sep="\t")

pairwise_q2_exonic <- data.table::fread("opt$counts_file", sep="\t", header = F, blank.lines.skip = T, 
  col.names = c("species1","species2","ancestor","transition_sp1","transition_sp2","Pattern","Count"))

pairwise_q2_exonic_ext_ages_sp1 <- pairwise_q2_exonic %>% 
    dplyr::select(species1, species2, ancestor, Pattern, transition_sp1, Count) %>% 
    mutate(joiner_colname = paste(species1, species2, ancestor, Pattern, transition_sp1, sep="__")) %>%
    group_by(joiner_colname) %>%
    summarise(Count_sp1 = sum(Count)) %>% as.data.frame() %>%
    dplyr::select(joiner_colname, Count_sp1)

pairwise_q2_exonic_ext_ages_sp2 <- pairwise_q2_exonic %>% 
    dplyr::select(species1, species2, ancestor, Pattern, transition_sp2, Count) %>% 
    mutate(joiner_colname = paste(species1, species2, ancestor, Pattern, transition_sp2, sep="__")) %>%
    group_by(joiner_colname) %>%
    summarise(Count_sp2 = sum(Count)) %>% as.data.frame() %>%
    dplyr::select(joiner_colname, Count_sp2)

pairwise_q2_exonic_ext_ages_comb <- pairwise_q2_exonic_ext_ages_sp1 %>%
    dplyr::left_join(pairwise_q2_exonic_ext_ages_sp2, 
                     by= "joiner_colname", suffix = c(".sp1", ".sp2")) %>% 
  tidyr::separate(joiner_colname, 
                      c("species1", "species2", "ancestor", "Pattern", "Transition"), sep = "__") %>%
  group_by(species1, species2, ancestor, Pattern) %>%
  mutate(Totsites_sp1 = sum(Count_sp1), Totsites_sp2 = sum(Count_sp2)) %>%
  as.data.frame() 

pairwise_q2_exonic_ext_ages_comb_det <- pairwise_q2_exonic_ext_ages_comb %>% 
  mutate(ref_5 = substring(Pattern, 1, 1), spc_5 = substring(Transition, 1, 1), 
         ref_base= substring(Pattern, 2, 2), spc_base = substring(Transition, 2, 2),
         ref_3 = substring(Pattern, 3, 3), spc_3 = substring(Transition, 3, 3)) %>%
  dplyr::left_join(seb_ages %>% dplyr::select(species, lifespan), by=c("species1"="species")) %>%
  dplyr::left_join(seb_ages %>% dplyr::select(species, lifespan), by=c("species2"="species"),
                   suffix =c("_sp1", "_sp2"))

pairwise_q2_exonic_ext_ages_comb1_classes <- pairwise_q2_exonic_ext_ages_comb_det %>% 
  inner_join(transition_table, by=c("ref_base"="base_orig", "spc_base"="base_transition")) %>%
  mutate(Mutation = paste0(ref_5, Class, ref_3)) %>%  
  group_by(species1, species2, ancestor, Mutation) %>%
  mutate(Class_sum_sp1 = sum(Count_sp1), Class_sum_sp2 = sum(Count_sp2)) %>%
  as.data.frame() %>% 
  dplyr::select(species1, species2, ancestor, Class, Mutation,
                Count_sp1, Count_sp2,
                Class_sum_sp1, Class_sum_sp2, lifespan_sp1, lifespan_sp2) %>%
  distinct()

pairwise_q2_exonic_ext_ages_comb1_classes_sum <- pairwise_q2_exonic_ext_ages_comb1_classes %>% 
  group_by(species1, species2, ancestor, Class) %>%
  mutate(prime_total_sp1 = sum(Class_sum_sp1), prime_total_sp2 = sum(Class_sum_sp2), 
         Class_prop_sp1 = Class_sum_sp1/prime_total_sp1, Class_prop_sp2 = Class_sum_sp2/prime_total_sp2) %>%
  as.data.frame() %>%
  group_by(species1, species2, ancestor) %>%
  mutate(prime_total_all_sp1=sum(prime_total_sp1), prime_total_all_sp2=sum(prime_total_sp2),
         Class_prop_all_sp1 = Class_sum_sp1/prime_total_all_sp1, 
         Class_prop_all_sp2 = Class_sum_sp2/prime_total_all_sp2,
         ratio_class_sum = Class_sum_sp1/Class_sum_sp2,
         ratio_class_prop = Class_prop_all_sp1/Class_prop_all_sp2,
         ratio_lifespan = lifespan_sp1/lifespan_sp2) %>%
  as.data.frame()

write.table(pairwise_q2_exonic_ext_ages_comb1_classes_sum, file = opt$totalsums_out, sep="\t", quote = F, row.names = F, col.names=T)

# Running the lm model
linear_reg_pairwise_q2_exonic <- pairwise_q2_exonic_ext_ages_comb1_classes_sum %>% 
  nest(data = c(-Mutation, -Class)) %>% 
  mutate(fit = map(data, ~ lm(ratio_class_sum ~ ratio_lifespan, data = .x)), 
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>% 
  dplyr::select(-fit, -data) %>% as.data.frame() %>% 
  arrange(Class) %>% filter(term != "(Intercept)")

write.table(linear_reg_pairwise_q2_exonic, file = opt$lm_outfile, sep="\t", quote = F, row.names = F, col.names=T)


