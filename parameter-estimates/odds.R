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
transition_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  as_log = FALSE
)

# Case 2: Gaining a function given a sibling did
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)

# Case 3: Speciation event
conditional_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1),
  i = 0, j = 0
)

# Case 4: Probability of preserving the function given sib lost
conditional_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = TRUE,
  array       = matrix(c(1,0), nrow = 1),
  i = 0, j = 0
)

# Case 5: Probability either gains a function
# This is equal to 1 - P(no gain at all)
1 - transition_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)


1 - transition_prob(
  p           = model2fit,
  duplication = FALSE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(0,0), nrow = 1),
  as_log = FALSE
)



