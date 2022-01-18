.PHONY: mcmc joint analysis-mcmc
mcmc-joint-geese: parameter-estimates/mcmc-joint-geese.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-geese.R \
		parameter-estimates/mcmc-joint-geese.Rout &

mcmc-joint-aphylo: parameter-estimates/mcmc-joint-aphylo.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-aphylo.R \
		parameter-estimates/mcmc-joint-aphylo.Rout &

# analysis-mcmc:
# 	Rscript --vanilla -e \
# 		'source("fig/mcmc-analysis.R", echo = TRUE)'
# analysis-mcmc-unif:
# 	Rscript --vanilla -e \
# 		'source("fig/mcmc-analysis-uniform-prior.R", echo = TRUE)'

# analysis-joint:
# 	Rscript --vanilla -e \
# 		'source("fig/mcmc-joint-analysis.R", echo = TRUE)'

.PHONY: ls
ls:
	ls -alt parameter-estimates/mcmc-PTH*.rds | wc -l ; \
		ls -alt parameter-estimates/mcmc-unif*.rds | wc -l
