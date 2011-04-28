## estimation functions should have
## INPUTS
##   ~ dat = data as a data.frame
##   ~ period.var = indicator of whether to include a period effect
##   ~ outcome.type = one of "gaussian", "binomial", "poisson"
##   ~ alpha = the type 1 error rate
## OUTPUTS
##   ~ a vector with three elements, in order:
##       [1] a point estimate for the treatment effect
##       [2] lower bound of (1-alpha) confidence interval
##       [3] lower bound of (1-alpha) confidence interval


fixed.effect <- function(dat, period.var, outcome.type, alpha) {
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}
	if(period.var==0){
		fit <- glm(y ~ trt + clust,
			   data=dat,
			   family=outcome.type,
			   offset=offsets)
	} else {
		fit <- glmer(y ~ trt + per + clust - 1,
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	}

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)$coef["trt",]
	ci <- est["Estimate"] + Z*est["Std. Error"]
	return(c(est["Estimate"], ci))

}

random.effect <- function(dat, period.var, outcome.type, alpha) {
	require(lme4)
	if(outcome.type=="poisson"){
		offsets <- log(dat[,"at.risk.time"])
	} else {
		offsets <- rep(0, nrow(dat))
	}

	if(period.var==0){
		fit <- glmer(y ~ trt + (1|clust),
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	} else {
		fit <- glmer(y ~ trt + per + (1|clust) - 1,
			     data=dat,
			     family=outcome.type,
			     offset=offsets)
	}

	## get CI
	Z <- qnorm(1-alpha/2)*c(-1,1)
	est <- summary(fit)@coefs["trt",]
	ci <- est[,"Estimate"] + Z*est[, "Std. Error"]
	return(c(est, ci))
}
