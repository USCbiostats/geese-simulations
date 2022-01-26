library(geese)

# Building the model
model2fit <- new_flock()
add_geese(
  p           = model2fit,
  annotations = list(0,0,0),
  geneid      = c(0,1,2),
  parent      = c(-1,0,0),
  duplication = c(TRUE, TRUE, TRUE)
  )

add_geese(
  p           = model2fit,
  annotations = list(0,0,0),
  geneid      = c(0,1,2),
  parent      = c(-1,0,0),
  duplication = !c(TRUE, TRUE, TRUE)
)

nfunctions <- 1L
# Building the model
term_overall_changes(model2fit, duplication = TRUE)  # Just constrain support
term_overall_changes(model2fit, duplication = FALSE) # Just constrain support

term_gains(model2fit, 0:(nfunctions - 1))
term_gains(model2fit, 0:(nfunctions - 1), FALSE)
term_loss(model2fit, 0:(nfunctions - 1))
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

# Indicator variable (this makes the difference)
term_k_genes_changing(model2fit, 1, TRUE)
term_k_genes_changing(model2fit, 1, FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

init_model(model2fit)

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

coefs <-
  c(
    `Overall changes at duplication` = 0,
    `Overall changes at speciation` = 0,
    `Only one gene changes at duplication` = -3.55887539187092,
    `Only one gene changes at speciation` = -2.45449172285053,
    `Gains 0 at duplication` = -0.433510761662447,
    `Gains 0 at speciation` = -4.10433924593962,
    `Loss 0 at duplication` = -1.78618476204973,
    `Loss 0 at speciation` = -5.68332000654832,
    `Root 1` = 6.29080043699408
  )

# Case 1: Moving from no function to both having a function

# Both gain
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  as_log = FALSE
)

# Only one gain
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(0,1), nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(0,1), nrow = 1),
  as_log = FALSE
)

# Either gain
1 - transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)

1 - transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)

# Both lose the function
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = TRUE,
  array       = matrix(c(1,1)*0, nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = TRUE,
  array       = matrix(c(1,1)*0, nrow = 1),
  as_log = FALSE
)

# Case 2: Gain given sib does not (so sib preserves a 0)

# Duplication
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  i = 0, j = 1
)

# Speciation
conditional_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = FALSE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 0
)

# Case 3: Preserve a one given sib does

# Duplication
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-9],
  state       = TRUE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)

# Speciation
conditional_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-9],
  state       = TRUE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)



