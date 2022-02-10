library(geese)
source("fig/plot_functions.R")

chains <- readRDS("parameter-estimates/mcmc-joint-geese-ver3.rds")
window(chains$mcmc[,-c(1,2)], start = end(chains$mcmc)*4/5) |>
  apply(2, quantile, probs = c(.025, .5, .975)) |>
  t()

# estimates <- window(chains$mcmc[,-c(1,2)], start = end(chains$mcmc)*4/5) |> colMeans()
estimates <- c(`Overall changes at duplication` = 0, `Overall changes at speciation` = 0,
               `Pairs of genes changing at duplication` = 1.07218091088764,
               `Pairs of genes changing at speciation` = 1.72931677452398, `Overall gains at duplication` = 2.24418138977199,
               `Overall gains at speciation` = 1.97323673341573, `Overall loses at duplication` = -3.32099148878171,
               `Overall loses at speciation` = 0.0300813296176934, `Root 1` = 0.0967000679510709,
               `Root 2` = 0.45371292244021)
# estimates[1:2] <- 0

traceplots(chains$mcmc[,-c(1,2)])

# Building the model
model2fit_odds <- new_flock()
add_geese(
  p           = model2fit_odds,
  annotations = list(c(0,0),c(0,0),c(0,0)),
  geneid      = c(0,1,2),
  parent      = c(-1,0,0),
  duplication = c(TRUE, TRUE, TRUE)
  )

add_geese(
  p           = model2fit_odds,
  annotations = list(c(0,0),c(0,0),c(0,0)),
  geneid      = c(0,1,2),
  parent      = c(-1,0,0),
  duplication = !c(TRUE, TRUE, TRUE)
)

term_overall_changes(model2fit_odds, duplication = TRUE)
term_overall_changes(model2fit_odds, duplication = FALSE)

term_pairwise_overall_change(model2fit_odds, duplication = 1)
term_pairwise_overall_change(model2fit_odds, duplication = 0)

term_overall_gains(model2fit_odds, 1)
term_overall_gains(model2fit_odds, 0)

term_overall_loss(model2fit_odds, 1)
term_overall_loss(model2fit_odds, 0)

rule_limit_changes(model2fit_odds, 0, 0, 4, TRUE)
rule_limit_changes(model2fit_odds, 1, 0, 4, FALSE)


init_model(model2fit_odds)
ans <- get_support(model2fit_odds)
ans <- do.call(rbind, ans)
View(ans)

# dat_geese_prior <- readRDS("parameter-estimates/mcmc-joint-geese.rds")
# colMeans(window(dat_geese_prior$mcmc, start = 15000))
# coefs <- c(
#   `Overall changes at duplication` = 0,
#   `Overall changes at speciation` = 0,
#   `Only one gene changes at duplication` = 4.98707582597087,
#   `Only one gene changes at speciation` = -3.03753881733255,
#   `Gains 0 at duplication` = -5.6549333970612,
#   `Gains 0 at speciation` = -5.23402105462181,
#   `Loss 0 at duplication` = -1.4447071742558,
#   `Loss 0 at speciation` = -6.71740250261233,
#   `Root 1` = 6.7549177338055
# )

# coefs <- c(
#   `Overall changes at duplication` = 0,
#   `Overall changes at speciation` = 0,
#   `Only one gene changes at duplication` = 4.00579512847272,
#   `Only one gene changes at speciation` = -3.18849622326337,
#   `Gains 0 at duplication` = -4.80443097143685,
#   `Gains 0 at speciation` = -5.32023091413009,
#   `Loss 0 at duplication` = -1.48689145790562,
#   `Loss 0 at speciation` = -6.89791021731092,
#   `Root 1` = 6.41661903218071
# )

coefs <- estimates # c(0,0, estimates)
# Case 1: Moving from no function to both having a function

# Both gain
transition_prob( # Duplication
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(rep(c(1,0), nfunctions), byrow = FALSE, nrow = nfunctions),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(rep(c(1,0), nfunctions), byrow = FALSE, nrow = nfunctions),
  as_log = FALSE
)

# Only one gain
transition_prob( # Duplication
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(c(1,rep(0, nfunctions * 2 - 1)), nrow = nfunctions),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(c(1,rep(0, nfunctions * 2 - 1)), nrow = nfunctions),
  as_log = FALSE
)

# Either gain
1 - transition_prob( # Duplication
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(rep(0, nfunctions * 2), nrow = nfunctions),
  as_log = FALSE
)

1 - transition_prob( # Speciation
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(FALSE, nfunctions),
  array       = matrix(rep(0, nfunctions * 2), nrow = nfunctions),
  as_log = FALSE
)

# Both lose the function
transition_prob( # Duplication
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(TRUE, nfunctions),
  array       = matrix(rep(0, nfunctions * 2), nrow = nfunctions),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-((length(estimates)-(nfunctions-1)):length(estimates))],
  state       = rep(TRUE, nfunctions),
  array       = matrix(rep(0, nfunctions * 2), nrow = nfunctions),
  as_log = FALSE
)

# Case 2: Gain given sib does not (so sib preserves a 0)

# Duplication
conditional_prob(
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 1
)

# Speciation
conditional_prob(
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 0
)

# Case 3: Preserve a zero given sib gained

# Duplication
conditional_prob(
  p           = model2fit_odds,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)

# Speciation
conditional_prob(
  p           = model2fit_odds,
  duplication = FALSE,
  params      = coefs[-9],
  state       = TRUE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)



