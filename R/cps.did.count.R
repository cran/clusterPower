#' Power simulations for cluster-randomized trials: Difference in Difference, Count Outcome.
#'
#' @description 
#' \loadmathjax
#'
#' This function utilizes iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs power simulations for difference in difference cluster randomized control trials using count outcomes
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per arm, between-cluster variance, 
#' two of the following: expected count in arm 1, expected count 
#' in arm 2, difference in counts between groups; significance level, 
#' analytic method, and whether or not progress updates should be displayed 
#' while the function is running.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters per arm; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified:
#' @param c1t0 Required. Expected outcome count in arm 1 at baseline.
#' Default is 0.
#' @param c2t0 Optional. Expected outcome count in arm 2 at baseline. If 
#' no quantity is provided, c2t0 = c1t0 is assumed.
#' @param c1t1 Optional. Expected outcome count in arm 1 at follow-up. 
#' If no quantity is provided, c1t1 = c1t0 is assumed.
#' @param c2t1 Required. Expected outcome count in arm 2 at follow-up.
#' @param c.diff Optional if c1t1 and c2t0 are provided. Expected difference 
#' in outcome count between groups, defined as 
#' c.diff = (c1t1 - c1t0) - (c2t1 - c2t0).
#' @param sigma_b_sq0 Pre-treatment (time == 0) between-cluster variance; 
#' accepts numeric scalar (indicating equal 
#' between-cluster variances for both arm) or a vector of length 2 specifying 
#' treatment-specific 
#' between-cluster variances
#' @param sigma_b_sq1 Post-treatment (time == 1) between-cluster variance; 
#' accepts numeric scalar (indicating equal 
#' between-cluster variances for both arm) or a vector of length 2 specifying 
#' treatment-specific 
#' between-cluster variances. For data simulation, sigma_b_sq1 is added to 
#' sigma_b_sq0, such that if sigma_b_sq0 = 5 
#' and sigma_b_sq1 = 2, the between-cluster variance at time == 1 equals 7. 
#' Default = 0.
#' @param alpha Significance level for power estimation, accepts value between 
#' 0 - 1; default = 0.05
#' @param family Distribution from which responses are simulated. Accepts Poisson 
#' ('poisson') or negative binomial ('neg.binom') (required); default = 'poisson'
#' @param analysis Family used for regression; currently only applicable for GLMM. 
#' Accepts c('poisson', 'neg.binom') (required); default = 'poisson'
#' @param negBinomSize Only used when generating simulated data from the 
#' negative binomial (family = 'neg.binom'), this is the target for number of 
#' successful trials, or the dispersion parameter (the shape parameter of the gamma 
#' mixing distribution). Must be strictly positive but need not be integer. 
#' Defaults to 1.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model 
#' (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') 
#' (required); default = 'glmm'
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time. Default = FALSE.
#' @param allSimData Option to output list of all simulated datasets. 
#' Default = FALSE
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge; default = FALSE.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' @param nofit Option to skip model fitting and analysis and only return the 
#' simulated data.
#' Default = \code{FALSE}.
#' @param seed Option to set the seed. Default is NA.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item Character string indicating total number of simulations, 
#'   distribution of simulated data, and regression family
#'   \item Number of simulations
#'   \item Data frame with columns 'Power' (Estimated statistical power), 
#'                'lower.95.ci' (Lower 95% confidence interval bound), 
#'                'upper.95.ci' (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Data frame containing families for distribution and analysis of simulated data
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting between-cluster variances at each time point for 
#'   each arm
#'   \item Vector containing expected counts and risk ratios based on user inputs
#'   \item Data frame with columns: 
#'                   'Period' (Pre/Post-treatment indicator), 
#'                   'Arm.2' (Arm indicator), 
#'                   'Value' (Mean response value)
#'   \item Data frame with columns: 
#'                   'Estimate' (Estimate of treatment effect for a given simulation), 
#'                   'Std.Err' (Standard error for treatment effect estimate), 
#'                   'Test.statistic' (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   'p.value', 
#'                   'converge' (Did simulated model converge?), 
#'                   'sig.val' (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
#'                   'y' (Simulated response value), 
#'                   'trt' (Indicator for arm), 
#'                   'clust' (Indicator for cluster), 
#'                   'period' (Indicator for time point)
#' }
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' 
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#' 
#' @examples 
#' 
#' # Estimate power for a trial with 7 clusters in both arms, those clusters having
#' # 9 subjects each, with sigma_b_sq0 = 0.1 in the first arm and 0.5 in the second arm. 
#' # We have estimated arm counts of 5 and 3 in the first and second arms, respectively, 
#' # and we use 100 simulated data sets analyzed by the GLMM method. The resulting 
#' # estimated power (if you set seed = 123) should be 0.86.
#' 
#' \dontrun{
#' did.count.sim = cps.did.count(nsim = 100, nsubjects = 9, nclusters = 7, 
#'                               c1t0 = 5, c1t1 = 5, c2t0 = 5, c2t1 = 8,  
#'                               sigma_b_sq0 = c(1, 0.5), sigma_b_sq1 = c(0.5, 0.8), 
#'                               family = 'poisson', analysis = 'poisson', 
#'                               method = 'glmm', seed = 123)
#' }
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#' @author Alexander R. Bogdan 
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @references Snjiders, T. & Bosker, R. Multilevel Analysis: an Introduction to Basic and Advanced Multilevel Modelling. London, 1999: Sage.
#' @references Elridge, S., Ukoumunne, O. & Carlin, J. The Intra-Cluster Correlation Coefficient in Cluster Randomized Trials: 
#' A Review of Definitions. International Statistical Review (2009), 77, 3, 378-394. doi: 10.1111/j.1751-5823.2009.00092.x
#' 
#'
#' @export

# Define function


cps.did.count = function(nsim = NULL,
                         nsubjects = NULL,
                         nclusters = NULL,
                         c1t0 = 0,
                         c2t0 = NULL,
                         c1t1 = NULL,
                         c2t1 = NULL,
                         c.diff = NULL,
                         sigma_b_sq0 = NULL,
                         sigma_b_sq1 = 0,
                         family = 'poisson',
                         analysis = 'poisson',
                         negBinomSize = 1,
                         method = 'glmm',
                         alpha = 0.05,
                         quiet = FALSE,
                         allSimData = FALSE,
                         poorFitOverride = FALSE,
                         lowPowerOverride = FALSE, 
                         timelimitOverride = TRUE,
                         seed = NA,
                         nofit = FALSE) {
  if (!is.na(seed)) {
    set.seed(seed = seed)
  }
  
  # Create vectors to collect iteration-specific values
  est.vector = vector("numeric", length = nsim)
  se.vector = vector("numeric", length = nsim)
  stat.vector = vector("numeric", length = nsim)
  pval.vector = vector("numeric", length = nsim)
  converge = vector("logical", length = nsim)
  values.vector = cbind(c(0, 0, 0, 0))
  simulated.datasets = list()
  
  # Create progress bar
  prog.bar =  progress::progress_bar$new(
    format = "(:spin) [:bar] :percent eta :eta",
    total = nsim,
    clear = FALSE,
    width = 100
  )
  prog.bar$tick(0)
  
  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5)
    abs(x - round(x)) < tol
  
  # Validate NSIM, NSUBJECTS, NCLUSTERS
  sim.data.arg.list = list(nsim, nclusters, nsubjects)
  sim.data.args = unlist(lapply(sim.data.arg.list, is.null))
  if (sum(sim.data.args) > 0) {
    stop("NSIM, NSUBJECTS, NCLUSTERS must all be specified. Please review your input values.")
  }
  min1.warning = " must be an integer greater than or equal to 1"
  if (!is.wholenumber(nsim) || nsim < 1) {
    stop(paste0("NSIM", min1.warning))
  }
  if (!is.wholenumber(nclusters) || nclusters < 1) {
    stop(paste0("NCLUSTERS", min1.warning))
  }
  if (!is.wholenumber(nsubjects) || nsubjects < 1) {
    stop(paste0("NSUBJECTS", min1.warning))
  }
  if (length(nclusters) > 2) {
    stop(
      "NCLUSTERS can only be a vector of length 1 (equal # of clusters per group) or 2 (unequal # of clusters per group)"
    )
  }
  # Set cluster sizes for arm 2 (if not already specified)
  if (length(nclusters) == 1) {
    nclusters[2] = nclusters[1]
  }
  # Set sample sizes for each cluster (if not already specified)
  if (length(nsubjects) == 1) {
    nsubjects[1:sum(nclusters)] = nsubjects
  }
  if (nclusters[1] == nclusters[2] &&
      length(nsubjects) == nclusters[1]) {
    nsubjects = rep(nsubjects, 2)
  }
  if (length(nsubjects) == 2 &&
      length(nclusters) == 2) {
    nsubjects <- c(rep(nsubjects[1], times = nclusters[1]),
      rep(nsubjects[2], times = nclusters[2]))
  }

  if (length(nclusters) == 2 &&
      length(nsubjects) != 1 &&
      length(nsubjects) != length(nclusters) &&
      length(nsubjects) != sum(nclusters)) {
    stop(
      "A cluster size must be specified for each cluster. If all cluster sizes are equal, please provide a single value for NSUBJECTS"
    )
  }
  
  # Validate C1, C2, C.DIFF
  parm1.arg.list = list(c1t0, c2t1, c.diff)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if (sum(parm1.args) > 1) {
    stop("At least two of the following terms must be specified: c1t0, c2t1, C.DIFF")
  }
  if (sum(parm1.args) == 0 && c.diff != abs(c2t1 - c1t0)) {
    stop("At least one of the following terms has been misspecified: c1t0, c2t1, C.DIFF")
  }
  # Set C1, C2, C.DIFF (if not already specified)
  
  if (!is.null(c.diff)) {
    if (is.null(c1t0)) {
      c1t0 = abs(c.diff - c2t1)
    }
    if (is.null(c2t1)) {
      c2t1 = abs(c1t0 - c.diff)
    }
  }
    if (is.null(c1t1)) {
      c1t1 = c1t0
    }
    if (is.null(c2t0)) {
      c2t0 = c1t0
    }
  
  if (is.null(c.diff)) {
    c.diff = (c1t1 - c1t0) - (c2t1 - c2t0)
  }
  
  # Validate sigma_b_sq0 & sigma_b_sq1
  sigma_b_sq.warning = " must be a scalar (equal between-cluster variance for both arms) or a vector of length 2,
  specifying between-cluster variances for each arm"
  if (!is.numeric(sigma_b_sq0) || any(sigma_b_sq0 < 0)) {
    stop("All values supplied to sigma_b_sq0 must be numeric values > 0")
  }
  if (!length(sigma_b_sq0) %in% c(1, 2)) {
    stop("sigma_b_sq0", sigma_b_sq.warning)
  }
  if (!length(sigma_b_sq1) %in% c(1, 2)) {
    stop("sigma_b_sq1", sigma_b_sq.warning)
  }
  if (!is.numeric(sigma_b_sq1) || any(sigma_b_sq1 < 0)) {
    stop("All values supplied to sigma_b_sq1 must be numeric values >= 0")
  }
  # Set sigma_b_sq0 & sigma_b_sq1 (if not already specified)
  if (length(sigma_b_sq0) == 1) {
    sigma_b_sq0[2] = sigma_b_sq0
  }
  if (length(sigma_b_sq1) == 1) {
    sigma_b_sq1[2] = sigma_b_sq1
  }
  sigma_b_sq1 = sigma_b_sq1 + sigma_b_sq0
  
  # Validate FAMILY, ANALYSIS, ALPHA, METHOD, QUIET, allSimData
  if (!is.element(family, c('poisson', 'neg.binom'))) {
    stop(
      "FAMILY must be either 'poisson' (Poisson distribution)
         or 'neg.binom'(Negative binomial distribution)"
    )
  }
  if (!is.element(analysis, c('poisson', 'neg.binom'))) {
    stop(
      "ANALYSIS must be either 'poisson' (Poisson regression)
         or 'neg.binom'(Negative binomial regression)"
    )
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  if (!is.element(method, c('glmm', 'gee'))) {
    stop(
      "METHOD must be either 'glmm' (Generalized Linear Mixed Model)
         or 'gee'(Generalized Estimating Equation)"
    )
  }
  if (!is.logical(quiet)) {
    stop(
      "QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)"
    )
  }
  if (!is.logical(allSimData)) {
    stop(
      "allSimData must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output"
    )
  }
  
  # Create indicators for PERIOD, TRT & CLUSTER
  period = rep(0:1, each = sum(nsubjects))
  trt = c(rep(1, length.out = sum(nsubjects[1:nclusters[1]])),
          rep(2, length.out = sum(nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])))
  clust = unlist(lapply(1:sum(nclusters), function(x)
    rep(x, length.out = nsubjects[x])))
  
  start.time = Sys.time()
  
  # Create simulation loop
  for (i in 1:nsim) {
    ## TIME == 0
    # Generate between-cluster effects for arm 1 and arm 2
    randint.ntrt.0 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq0[1]))
    randint.trt.0 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq0[2]))
    
    # Create arm 1 y-value
    y0.ntrt.intercept = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.ntrt.0[x], length.out = nsubjects[x])))
    y0.ntrt.linpred = y0.ntrt.intercept + log(c1t0)
    y0.ntrt.prob = exp(y0.ntrt.linpred)
    if (family == 'poisson') {
      y0.ntrt = stats::rpois(length(y0.ntrt.prob), y0.ntrt.prob)
    }
    if (family == 'neg.binom') {
      y0.ntrt = stats::rnbinom(length(y0.ntrt.prob), size = negBinomSize, mu = y0.ntrt.prob)
    }
    
    # Create arm 2 y-value
    y0.trt.intercept = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.trt.0[x], length.out = nsubjects[nclusters[1] + x])))
    y0.trt.linpred = y0.trt.intercept + log(c2t0)
    y0.trt.prob = exp(y0.trt.linpred)
    if (family == 'poisson') {
      y0.trt = stats::rpois(length(y0.trt.prob), y0.trt.prob)
    }
    if (family == 'neg.binom') {
      y0.trt = stats::rnbinom(length(y0.trt.prob), size = 1, mu = y0.trt.prob)
    }
    
    ## TIME == 1
    # Generate between-cluster effects for arm 1 and arm 2
    randint.ntrt.1 = stats::rnorm(nclusters[1], mean = 0, sd = sqrt(sigma_b_sq1[1]))
    randint.trt.1 = stats::rnorm(nclusters[2], mean = 0, sd = sqrt(sigma_b_sq1[2]))
    
    # Create arm 1 y-value
    y1.ntrt.intercept = unlist(lapply(1:nclusters[1], function(x)
      rep(randint.ntrt.1[x], length.out = nsubjects[x])))
    y1.ntrt.linpred = y1.ntrt.intercept + log(c1t1)
    y1.ntrt.prob = exp(y1.ntrt.linpred)
    if (family == 'poisson') {
      y1.ntrt = stats::rpois(length(y1.ntrt.prob), y1.ntrt.prob)
    }
    if (family == 'neg.binom') {
      y1.ntrt = stats::rnbinom(length(y1.ntrt.prob), size = 1, mu = y1.ntrt.prob)
    }
    
    # Create arm 2 y-value
    y1.trt.intercept = unlist(lapply(1:nclusters[2], function(x)
      rep(randint.trt.1[x], length.out = nsubjects[nclusters[1] + x])))
    y1.trt.linpred = y1.trt.intercept + log(c2t1)
    y1.trt.prob = exp(y1.trt.linpred)
    if (family == 'poisson') {
      y1.trt = stats::rpois(length(y1.trt.prob), y1.trt.prob)
    }
    if (family == 'neg.binom') {
      y1.trt = stats::rnbinom(length(y1.trt.prob), size = 1, mu = y1.trt.prob)
    }
    
    # Create single response vector
    y = c(y0.ntrt, y0.trt, y1.ntrt, y1.trt)
    
    # Create and store data for simulated dataset
    sim.dat = data.frame(
      y = y,
      trt = trt,
      period = period,
      clust = clust
    )
    if (allSimData == TRUE) {
      simulated.datasets = append(simulated.datasets, list(sim.dat))
    }
    
    # option to return simulated data only
    if (nofit == TRUE) {
      if (!exists("nofitop")) {
        nofitop <- data.frame(
          period = sim.dat['period'],
          cluster = sim.dat['clust'],
          arm = sim.dat['trt'],
          y1 = sim.dat["y"]
        )
      } else {
        nofitop[, length(nofitop) + 1] <- sim.dat["y"]
      }
      if (length(nofitop) == (nsim + 3)) {
        temp1 <- seq(1:nsim)
        temp2 <- paste0("y", temp1)
        colnames(nofitop) <- c('period', 'cluster', 'arm', temp2)
      }
      if (length(nofitop) != (nsim + 3)) {
        next()
      }
      return(nofitop)
    }
    
    
    # Calculate mean values for given simulation
    iter.values = cbind(stats::aggregate(y ~ trt + period, data = sim.dat, mean)[, 3])
    values.vector = values.vector + iter.values
    
    # Set warnings to OFF
    # Note: Warnings will still be stored in 'warning.list'
    options(warn = -1)
    
    # Fit GLMM (lmer)
    if (method == 'glmm') {
      if (analysis == 'poisson') {
        my.mod = lme4::glmer(
          y ~ trt + period + trt:period + (1 |
                                             clust),
          data = sim.dat,
          family = stats::poisson(link = 'log')
        )
      }
      if (analysis == 'neg.binom') {
        my.mod = lme4::glmer.nb(y ~ trt + period + trt:period + (1 |
                                                                   clust), data = sim.dat)
      }
      glmm.values = summary(my.mod)$coefficient
      est.vector[i] = glmm.values['trt:period', 'Estimate']
      se.vector[i] = glmm.values['trt:period', 'Std. Error']
      stat.vector[i] = glmm.values['trt:period', 'z value']
      pval.vector[i] = glmm.values['trt:period', 'Pr(>|z|)']
      converge[i] = is.null(my.mod@optinfo$conv$lme4$messages)
    }
    
    # Set warnings to ON
    options(warn = 0)
    
    # Fit GEE (geeglm)
    if (method == 'gee') {
      sim.dat = dplyr::arrange(sim.dat, clust)
      my.mod = geepack::geeglm(
        y ~ trt + period + trt:period,
        data = sim.dat,
        family = stats::poisson(link = 'log'),
        id = clust,
        corstr = "exchangeable"
      )
      gee.values = summary(my.mod)$coefficients
      est.vector[i] = gee.values['trt:period', 'Estimate']
      se.vector[i] = gee.values['trt:period', 'Std.err']
      stat.vector[i] = gee.values['trt:period', 'Wald']
      pval.vector[i] = gee.values['trt:period', 'Pr(>|W|)']
      converge[i] <- ifelse(summary(my.mod)$error == 0, TRUE, FALSE)
    }
    
    # option to stop the function early if fits are singular
    if (poorFitOverride == FALSE && converge[i] == FALSE) {
      if (sum(converge == FALSE, na.rm = TRUE) > (nsim * .25)) {
        stop(
          "more than 25% of simulations are singular fit: check model specifications"
        )
      }
    }
    
    # stop the loop if power is <0.5
    if (lowPowerOverride == FALSE && i > 50 && (i %% 10 == 0)) {
      sig.val.temp <-
        ifelse(pval.vector < alpha, 1, 0)
      pval.power.temp <- sum(sig.val.temp, na.rm = TRUE) / i
      if (pval.power.temp < 0.5) {
        stop(
          paste(
            "Calculated power is < ",
            pval.power.temp,
            ". Set lowPowerOverride == TRUE to run the simulations anyway.",
            sep = ""
          )
        )
      }
    }
    
    # Update progress information
    if (quiet == FALSE) {
      if (i == 1) {
        avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
        time.est = avg.iter.time * (nsim - 1) / 60
        hr.est = time.est %/% 60
        min.est = round(time.est %% 60, 3)
        if (min.est > 2 && timelimitOverride == FALSE){
          stop(paste0("Estimated completion time: ",
                      hr.est,
                      'Hr:',
                      min.est,
                      'Min'
          ))
        }
        message(
          paste0(
            'Begin simulations :: Start Time: ',
            Sys.time(),
            ' :: Estimated completion time: ',
            hr.est,
            'Hr:',
            min.est,
            'Min'
          )
        )
      }
      # Iterate progress bar
      prog.bar$update(i / nsim)
      Sys.sleep(1 / 100)
      
      if (i == nsim) {
        message(paste0("Simulations Complete! Time Completed: ", Sys.time()))
      }
    }
  }
  
  ## Output objects
  # Create object containing summary statement
  summary.message = paste0(
    "Monte Carlo Power Estimation based on ",
    nsim,
    " Simulations: Difference in Difference Design, Count Outcome\nData Simulated from ",
    switch(family, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
    " distribution\nAnalyzed using ",
    switch(analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'),
    " regression"
  )
  
  # Create method object
  long.method = switch(method, glmm = 'Generalized Linear Mixed Model',
                       gee = 'Generalized Estimating Equation')
  
  # Store simulation output in data frame
  cps.model.est = data.frame(
    Estimate = as.vector(unlist(est.vector)),
    Std.err = as.vector(unlist(se.vector)),
    Test.statistic = as.vector(unlist(stat.vector)),
    p.value = as.vector(unlist(pval.vector)),
    converge = converge
  )
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  cps.model.temp <- dplyr::filter(cps.model.est, converge == TRUE)
  power.parms <- confintCalc(alpha = alpha,
                             nsim = nsim,
                             p.val = cps.model.temp[, 'p.value'])
  
  # Create object containing inputs
  c1.c2.rr = round(exp(log(c1t1) - log(c2t1)), 3)
  c2.c1.rr = round(exp(log(c2t1) - log(c1t1)), 3)
  inputs = t(data.frame(
    'Arm.1' = c("count" = abs(c1t1), "risk.ratio" = c1.c2.rr),
    'Arm.2' = c("count" = abs(c2t1), 'risk.ratio' = c2.c1.rr),
    'Difference' = c("count" = c.diff, 'risk.ratio' = c2.c1.rr - c1.c2.rr)
  ))
  
  # Create object containing arm & time-specific differences
  values.vector = values.vector / nsim
  differences = data.frame(
    Period = c(0, 0, 1, 1),
    Arm.2 = c(0, 1, 0, 1),
    Values = round(values.vector, 3)
  )
  
  # Create object containing group-specific cluster sizes
  cluster.sizes = list('Arm.1' = nsubjects[1:nclusters[1]],
                       'Arm.2' = nsubjects[(nclusters[1] + 1):(nclusters[1] + nclusters[2])])
  
  # Create object containing number of clusters
  n.clusters = t(data.frame(
    "Arm.1" = c("n.clust" = nclusters[1]),
    "Arm.2" = c("n.clust" = nclusters[2])
  ))
  
  # Create object containing group-specific variance parameters
  var.parms = list(
    "Time.Point.0" = data.frame(
      'Arm.1' = c("sigma_b_sq" = sigma_b_sq0[1]),
      'Arm.2' = c("sigma_b_sq" = sigma_b_sq0[2])
    ),
    "Time.Point.1" = data.frame(
      'Arm.1' = c("sigma_b_sq" = sigma_b_sq1[1]),
      'Arm.2' = c("sigma_b_sq" = sigma_b_sq1[2])
    )
  )
  
  # Create object containing FAMILY & REGRESSION parameters
  dist.parms = rbind(
    'Family:' = paste0(switch(
      family, poisson = 'Poisson', neg.binom = 'Negative Binomial'
    ), ' distribution'),
    'Analysis:' = paste0(switch(
      analysis, poisson = 'Poisson', neg.binom = 'Negative Binomial'
    ), ' distribution')
  )
  colnames(dist.parms) = "Data Simuation & Analysis Parameters"
  
  # Create list containing all output and return
  complete.output = structure(
    list(
      "overview" = summary.message,
      "nsim" = nsim,
      "power" = power.parms,
      "method" = long.method,
      "alpha" = alpha,
      "cluster.sizes" = cluster.sizes,
      "n.clusters" = n.clusters,
      "variance.parms" = var.parms,
      "dist.parms" = dist.parms,
      "inputs" = inputs,
      "differences" = differences,
      "model.estimates" = cps.model.est,
      "sim.data" = simulated.datasets
    ),
    class = 'crtpwr'
  )
  
  return(complete.output)
}
