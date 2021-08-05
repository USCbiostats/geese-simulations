.PHONY: mcmc joint analysis-mcmc
mcmc: parameter-estimates/mcmc.R
	R CMD BATCH --vanilla parameter-estimates/mcmc.R \
		parameter-estimates/mcmc.Rout &
mcmc-no-prior: parameter-estimates/mcmc-no-prior.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-no-prior.R \
		parameter-estimates/mcmc-no-prior.Rout &


joint: parameter-estimates/mcmc-joint.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint.R \
		parameter-estimates/mcmc-joint.Rout &

secondrun: parameter-estimates/mcmc-joint-geese-second-run.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-geese-second-run.R \
		parameter-estimates/mcmc-joint-geese-second-run.Rout &

analysis-mcmc:
	Rscript --vanilla -e \
		'source("parameter-estimates/mcmc-analysis.R", echo = TRUE)'

analysis-joint:
	Rscript --vanilla -e \
		'source("parameter-estimates/mcmc-joint-analysis.R", echo = TRUE)'

