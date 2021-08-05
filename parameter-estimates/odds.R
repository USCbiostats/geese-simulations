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
term_overall_changes(model2fit, duplication = TRUE)
term_overall_changes(model2fit, duplication = FALSE)
term_genes_changing(model2fit, duplication = TRUE)
term_gains(model2fit, 0:(nfunctions - 1))
term_gains(model2fit, 0:(nfunctions - 1), FALSE)
term_loss(model2fit, 0:(nfunctions - 1), FALSE)

rule_limit_changes(model2fit, 0, 0, 4, TRUE)
rule_limit_changes(model2fit, 1, 0, 4, FALSE)

init_model(model2fit)

# dat_geese_prior <- readRDS("parameter-estimates/mcmc-joint-geese-prior.rds")
# colMeans(window(dat_geese_prior$mcmc, start = 1e4))
coefs <- c(
  `Overall changes at duplication` = 0,
  `Overall changes at speciation` = 0,
  `Num. of genes changing at duplication` = -1.80204158815901,
  `Gains 0 at duplication` = 0.0768957326737147,
  `Gains 0 at speciation` = -2.15269078801566,
  `Loss 0 at speciation` = -3.42075856213656,
  `Root 1` = 1.64267185610882
)

# Case 1: Moving from no function to both having a function

# Both gain
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  as_log = FALSE
)

# Only one gain
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,1), nrow = 1),
  as_log = FALSE
)

# Either gain
1 - transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)

1 - transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)

# Both lose the function
transition_prob( # Duplication
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = TRUE,
  array       = matrix(c(1,1)*0, nrow = 1),
  as_log = FALSE
)

transition_prob( # Speciation
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = TRUE,
  array       = matrix(c(1,1)*0, nrow = 1),
  as_log = FALSE
)

# Case 2: Gain given sib does not (so sib preserves a 0)

# Duplication
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 0
) * 2

# Speciation
conditional_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 0
)

# Case 3: Preserve a one given sib does

# Duplication
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = TRUE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)

# Speciation
conditional_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = TRUE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)



