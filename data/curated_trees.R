library(aphylo)
library(geese)

list.files("data-raw/pthr16LimitedCandTreeSet")


trees <- readRDS("data-raw/pthr16LimitedCandTreeSet/candtreesanno>1pos+neg.rds")

model1 <- aphylo_to_geese(trees[[1]])
term_overall_changes(model1, event_type = 2)
term_k_genes_changing(model1, 1, event_type = 2)
term_gains(model1, 0:1, event_type = 0)
term_gains(model1, 0:1, event_type = 1)
term_loss(model1, 0:1, event_type = 0)
term_loss(model1, 0:1, event_type = 1)

rule_limit_changes(model1, 0, lb = 0, ub = 4, event_type = 2)
init_model(model1)


geese_mle(model1)
