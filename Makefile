.PHONY: mcmc joint analysis-mcmc
mcmc: parameter-estimates/mcmc.R
	R CMD BATCH --vanilla parameter-estimates/mcmc.R \
		parameter-estimates/mcmc.Rout &
mcmc-unif-prior: parameter-estimates/mcmc-unif-prior.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-unif-prior.R \
		parameter-estimates/mcmc-unif-prior.Rout &


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

