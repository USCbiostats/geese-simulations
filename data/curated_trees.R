library(aphylo)
library(geese)

list.files("data-raw/pthr16LimitedCandTreeSet")


trees <- readRDS("data-raw/pthr16LimitedCandTreeSet/candtreesanno>1pos+neg.rds")

for (i in 1:length(trees)) {

  model1 <- aphylo_to_geese(trees[[i]])
  parse_polytomies(model1)
  term_overall_changes(model1, event_type = 2)
  term_k_genes_changing(model1, 1, event_type = 2)
  term_gains(model1, 0:1, event_type = 0)
  term_gains(model1, 0:1, event_type = 1)
  term_loss(model1, 0:1, event_type = 0)
  term_loss(model1, 0:1, event_type = 1)

  if (parse_polytomies(model1) > 10)
    next

  rule_limit_changes(model1, 0, lb = 0, ub = 4, event_type = 2)
  stop()
}


init_model(model1)

ans_mle  <- geese_mle(model1)
ans_mcmc  <- geese_mcmc(
  model1, nsteps=2000,
  kernel = fmcmc::kernel_ram(
    lb = -10, ub = 10, fixed = c(TRUE, rep(FALSE, 11))
    )
  )
ans_pred <- predict_geese(model1, par = ans_mle$par, leave_one_out = TRUE, only_annotated = TRUE)
