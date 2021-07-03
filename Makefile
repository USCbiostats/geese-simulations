.PHONY: mcmc joint analysis-mcmc
mcmc: parameter-estimates/mcmc.R
	R CMD BATCH --vanilla parameter-estimates/mcmc.R \
		parameter-estimates/mcmc.Rout &
mcmc-no-prior: parameter-estimates/mcmc-no-prior.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-no-prior.R \
		parameter-estimates/mcmc-no-prior.Rout &


joint: parameter-estimates/mcmc-joint-geese.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-geese.R \
		parameter-estimates/mcmc-joint-geese.Rout &

analysis-mcmc:
	Rscript --vanilla -e \
		'source("parameter-estimates/mcmc-analysis.R", echo = TRUE)'
