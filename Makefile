.PHONY: mcmc joint analysis-mcmc
mcmc-joint-geese-ver2: parameter-estimates/mcmc-joint-geese-ver2.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-geese-ver2.R \
		parameter-estimates/mcmc-joint-geese-ver2.Rout &

mcmc-joint-aphylo: parameter-estimates/mcmc-joint-aphylo.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-joint-aphylo.R \
		parameter-estimates/mcmc-joint-aphylo.Rout &

mcmc-unif-prior-curated: parameter-estimates/mcmc-unif-prior-curated.R
	R CMD BATCH --vanilla parameter-estimates/mcmc-unif-prior-curated.R \
		parameter-estimates/mcmc-unif-prior-curated-v2.Rout &



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
