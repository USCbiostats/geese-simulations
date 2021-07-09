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

transition_prob(
  p           = model2fit,
  duplication = TRUE,
  params      = coefs[-7],
  state       = FALSE,
  array       = matrix(c(1,1), nrow = 1), as_log = FALSE
)

coefs <- coefs[-c(1:2, 7)]


# Cases: (Spec, Dupl) x (0, 1) x (Gain, Lose) x (Sib 0, Sib 1) -----------------

# Case: (Spec, 0, Gain, Sib 0), Gain a function in spec given sib doesn't have
# - (A has function) - (No change)
# - c(0, 0, 1, 0) - 0
counts_1 <- (c(0, 0, 1, 0) - 0) * coefs
1/(1 + exp(-sum(counts_1)))

# Case: (Dupl, 0, Gain, Sib 0), Gain a function in dupl given sib doesn't have
# - (A has function) - (No change)
# - c(1, 1, 0, 0) - 0
counts_1 <- (c(1, 1, 0, 0) - 0) * coefs
1/(1 + exp(-sum(counts_1)))


# Case 1: Only B gains the function
# - (B has the function) - (Both have the function)
# - c(1, 1, 0, 0) - c(2, 2, 0, 0)
counts_1 <- (c(1, 1, 0, 0) - c(2, 2, 0, 0)) * coefs
1/(1 + exp(-sum(counts_1)))

# Case 2: No changes at spec event (given B)
# - (No changes) - (A lost a func)
counts_1 <- (c(0, 0, 0, 0) - c(0, 0, 0, 1) ) * coefs
1/(1 + exp(-sum(counts_1)))

# Case 3: A gains a function from scratch (given B doesn't)
# - (Neither has a function) - (A gains a function)
counts_1 <- (c(0,0,0,0) - c(0, 0, 1, 0)) * coefs
1/(1 + exp(-sum(counts_1)))

# Case 4: A preserves the function given B lost it
counts_1 <- c(-1, 0, 0, 0) * coefs
1 - 1/(1 + exp(-sum(counts_1)))

# Case 5: Likelihood of either gaining the function
# = P(A, ~B|~P) + P(~A, B|~P) + P(A, B|~P)
# = P(A|~B,~P)P(~B|P) + P(~B|A,~P)P(~A|P) +
