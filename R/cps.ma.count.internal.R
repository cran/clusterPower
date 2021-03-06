#' Model fits for simulations for multi-arm designs with count outcome.
#'
#' Generally called from \code{cps.ma.count()}, this function uses iterative
#' simulations to model significance of treatment effects for cluster-randomized
#' controlled trials. Users can modify a variety of parameters to suit the
#' simulations to their desired experimental situation.
#'
#' This function can be called directly in order to give the user access to the
#' simulated model fits in addition to the simulated data, the latter of which
#' can also be accessed here or using the function \code{cps.ma.count()}. This
#' function does not produce power estimates, just simulated data and model fits.
#'
#' To call this function independenty, users must specify the desired number of
#' simulations, number of subjects per cluster, number of clusters per treatment
#' arm, group proportions, and between-cluster variance. Significance level,
#' analytic method, progress updates and simulated data set output may also be
#' specified.
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#'
#' @param family A string consisting of either 'poisson' or 'neg.binom'.
#'
#' @param str.nsubjects Number of subjects per treatment group; accepts a list
#' with one entry per arm. Each entry is a vector containing the number of
#' subjects per cluster (required).
#'
#' @param counts Expected outcome for each arm; accepts a vector
#' of length \code{narms} (required).
#'
#' @param sigma_b_sq Between-cluster variance; accepts a vector of length
#' \code{narms} (required).
#'
#' @param alpha Significance level; default = 0.05.
#'
#' @param method Analytical method, either Generalized Linear Mixed Effects
#' Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm',
#' 'gee') (required); default = 'glmm'.
#'
#' @param analysis Family used for data analysis; currently only applicable when \code{method = 'glmm'}.
#' Accepts c('poisson', 'neg.binom'); default = 'poisson'. Required.
#'
#' @param negBinomSize Only used when generating simulated data from the
#' negative binomial (family = 'neg.binom'), this is the target for number of
#' successful trials, or the dispersion parameter (the shape parameter of the gamma
#' mixing distribution). Must be positive and defaults to 1. Required when
#' family = 'neg.binom'.
#'
#' @param quiet When set to FALSE, displays simulation progress and estimated
#' completion time; default is FALSE.
#'
#' @param all.sim.data Option to output list of all simulated datasets;
#' default = FALSE.
#'
#' @param seed Option to set.seed. Default is NULL.
#'
#' @param poor.fit.override Option to override \code{stop()} if more than 25\%
#' of fits fail to converge.
#'
#' @param low.power.override Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#'
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated
#' completion time is more than 2 minutes. Defaults to TRUE.
#'
#' @param tdist Logical; use t-distribution instead of normal distribution for
#' simulation values, default = FALSE.
#'
#' @param cores A string ("all") NA, or numeric value indicating the number of
#' cores to be used for parallel computing. When this option is set to NA, no
#' parallel computing is used.
#'
#' @param nofit Option to skip model fitting and analysis and return the
#' simulated data. Defaults to \code{FALSE}.
#'
#' @param opt Option to fit with a different optimizer algorithm. Setting this
#' to "auto" tests an example fit using the \code{nloptr} package and selects
#' the first algorithm that converges.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item List of length(nsim) containing gee- or glmm-fitted the model
#'   summaries.
#'   \item Compares fitted model to a model for H0 using ML (anova).
#'   \item List of data frames, each containing:
#'                   "y" (Simulated response value),
#'                   "trt" (Indicator for treatment group),
#'                   "clust" (Indicator for cluster)
#'   \item A vector of length \code{nsim} consisting of 1 and 0.
#'           When a model fails to converge, failed.to.converge==1, otherwise 0.
#' }
#'
#' @examples
#' \dontrun{
#'
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' counts.example <- c(75, 120, 100)
#' sigma_b_sq.example <- c(0.2, 0.1, 0.1)
#'
#' count.ma.rct <- cps.ma.count.internal (nsim = 10,
#'                             str.nsubjects = nsubjects.example,
#'                             counts = counts.example,
#'                             sigma_b_sq = sigma_b_sq.example,
#'                             alpha = 0.05, all.sim.data = FALSE,
#'                             seed = 123, cores="all", poor.fit.override=TRUE,
#'                             opt = "NLOPT_LN_BOBYQA")
#' }
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @noRd
cps.ma.count.internal <-
  function(nsim = 1000,
           str.nsubjects = NULL,
           counts = NULL,
           family = "poisson",
           analysis = "poisson",
           negBinomSize = 1,
           sigma_b_sq = NULL,
           alpha = 0.05,
           quiet = FALSE,
           method = 'glmm',
           all.sim.data = FALSE,
           seed = NA,
           poor.fit.override = FALSE,
           low.power.override = FALSE,
           timelimitOverride = TRUE,
           tdist = FALSE,
           cores = cores,
           nofit = FALSE,
           opt = opt) {
    # Create vectors to collect iteration-specific values
    simulated.datasets <- list()
    goodopt <- opt

    # Create NCLUSTERS, NARMS, from str.nsubjects
    narms = length(str.nsubjects)
    nclusters = sapply(str.nsubjects, length)
    
    # This container keeps track of how many models failed to converge
    converged <- rep(FALSE, nsim)
    
    # Create a container for the simulated.dataset and model output
    sim.dat = vector(mode = "list", length = nsim)
    model.values <- list()
    model.compare <- list()
    
    # option for reproducibility
    if (!is.na(seed)) {
      set.seed(seed = seed)
    }
    
    # Create indicators for treatment group & cluster for the sim.data output
    trt1 = list()
    clust1 = list()
    index <- 0
    for (arm in 1:length(str.nsubjects)) {
      trt1[[arm]] = list()
      clust1[[arm]] =  list()
      for (cluster in 1:length(str.nsubjects[[arm]])) {
        index <- index + 1
        trt1[[arm]][[cluster]] = rep(arm, sum(str.nsubjects[[arm]][[cluster]]))
        clust1[[arm]][[cluster]] = rep(index, sum(str.nsubjects[[arm]][[cluster]]))
      }
    }
    
    #Alert the user if using t-distribution
    if (tdist == TRUE) {
      print("using t-distribution because tdist = TRUE")
    }
    #make the simulated data
    trt <-  as.factor(unlist(trt1))
    clust <- as.factor(unlist(clust1))
    if (length(trt) != length(clust)) {
      stop("trt and clust are not the same length, see line 134")
    }
    gc()
    sim.dat <- matrix(nrow = length(clust), ncol = nsim)
    
    # function to produce the simulated data
    make.sim.dat <- function(tdist. = tdist,
                             counts. = counts,
                             nclusters. = nclusters,
                             sigma_b_sq. = sigma_b_sq,
                             str.nsubjects. = str.nsubjects,
                             family. = family) {
      # Generate between-cluster effects for non-treatment and treatment
      if (tdist. == TRUE) {
        randint = mapply(function(n, df)
          stats::rt(n, df = df),
          n = nclusters.,
          df = Inf)
      } else {
        randint = mapply(
          function(nc, s, mu)
            stats::rnorm(nc, mean = mu,
                         sd = sqrt(s)),
          nc = nclusters.,
          s = sigma_b_sq.,
          mu = 0
        )
      }
      if (typeof(randint) == "list") {
        randint.holder <- list()
        for (j in 1:length(counts.)) {
          randint.holder[[j]] <- log(counts.[j]) + randint[[j]]
        }
        randintrandint <- sapply(randint.holder, exp)
      } else {
        randint.holder <-
          matrix(nrow = nclusters.[1], ncol = length(counts.))
        for (j in 1:length(counts.)) {
          randint.holder[, j] <- log(counts.[j]) + randint[, j]
        }
        randintrandint <- exp(randint.holder)
      }
      
      # Create y-value
      y.intercept <- vector(mode = "numeric",
                            length = sum(unlist(str.nsubjects.)))
      y.intercept <- sapply(1:sum(nclusters.),
                            function(x)
                              rep(unlist(randintrandint)[x],
                                  length.out = unlist(str.nsubjects.)[x]))
      y.intercept <- as.vector(unlist(y.intercept))
      if (family. == 'poisson') {
        y <-  stats::rpois(length(y.intercept), y.intercept)
      }
      if (family. == 'neg.binom') {
        y <-
          stats::rnbinom(length(y.intercept), size = negBinomSize, mu = y.intercept)
      }
      return(y)
    }
    sim.dat <- replicate(
      nsim,
      make.sim.dat(
        tdist. = tdist,
        nclusters. = nclusters,
        sigma_b_sq. = sigma_b_sq,
        str.nsubjects. = str.nsubjects,
        counts. = counts,
        family. = family
      )
    )
    
    #option to return simulated data only
    if (nofit == TRUE) {
      sim.dat <- data.frame(trt, clust, sim.dat)
      sim.num <- 1:nsim
      temp <- paste0("y", sim.num)
      colnames(sim.dat) <- c("arm", "cluster", temp)
      return(sim.dat)
    }
    
    `%fun%` <- foreach::`%dopar%`
    if (is.na(cores)) {
      `%fun%` <- foreach::`%do%`
    }
    
    #setup for parallel computing
    ## Do computations with multiple processors:
    ## Number of cores:
    if (!is.na(cores)) {
      if (cores == "all") {
        nc <- parallel::detectCores()
      } else {
        nc <- cores
      }
      ## Create clusters and initialize the progress bar
      cl <-
        parallel::makeCluster(rep("localhost", nc), type = "SOCK")
      doParallel::registerDoParallel(cl)
    }
    pb <- txtProgressBar(min = 1, max = nsim, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # Update simulation progress information
    sim.start <- Sys.time()
    lme4::glmer(sim.dat[, 1] ~ as.factor(trt) + (1 |
                                                   clust),
                family = stats::poisson(link = 'log'))
    avg.iter.time = as.numeric(difftime(Sys.time(), sim.start, units = 'secs'))
    time.est = avg.iter.time * (nsim - 1) / 60
    hr.est = time.est %/% 60
    min.est = round(time.est %% 60, 3)
    #time limit override (for Shiny)
    if (min.est > 2 && timelimitOverride == FALSE) {
      stop(paste0("Estimated completion time: ",
                  hr.est,
                  'Hr:',
                  min.est,
                  'Min'))
    }
    if (quiet == FALSE) {
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
      # initialize progress bar
      if (is.na(cores)) {
        prog.bar =  progress::progress_bar$new(
          format = "(:spin) [:bar] :percent eta :eta",
          total = 5,
          clear = FALSE,
          show_after = 0
        )
        prog.bar$tick(0)
      }
    }
    
    if (is.na(cores) & quiet == FALSE) {
      # Iterate progress bar
      prog.bar$update(1 / 5)
      Sys.sleep(1 / 100)
    }
    
    if (method == "glmm") {
      # Fit the models
      
      if (!is.na(cores) & quiet == FALSE) {
        message("Fitting models")
      }
      
      if (analysis == 'poisson') {
        if (opt == "auto") {
          message("testing optimizer algorithms:")
          algoptions <- c(
            "NLOPT_LN_BOBYQA",
            "NLOPT_GN_CRS2_LM",
            "NLOPT_LN_COBYLA",
            "NLOPT_LN_NEWUOA",
            "NLOPT_LN_NEWUOA_BOUND",
            "NLOPT_LN_NELDERMEAD",
            "NLOPT_LN_SBPLX"
          )
          for (i in 1:length(algoptions)) {
            print(algoptions[i])
            R.utils::withTimeout(Sys.sleep(10), timeout = 30)
            model_flex1 <-
              try(lme4::glmer(
                sim.dat[, 1] ~ as.factor(trt) + (1 | clust),
                family = stats::poisson(link = 'log'),
                control = lme4::glmerControl(
                  optimizer = "nloptwrap",
                  optCtrl = list(
                    algorithm = algoptions[i],
                    maxeval = 1e7,
                    xtol_abs = 1e-9,
                    ftol_abs = 1e-9
                  )
                )
              ))
            if (class(model_flex1) != "try-error") {
              if (is.null(model_flex1@optinfo$conv$lme4$messages)) {
                goodopt <- algoptions[i]
                break
              }
            }
          }
        }
        my.mod <- foreach::foreach(
          i = 1:nsim,
          .options.parallel = opts,
          .packages = c("lme4", "optimx", "nloptr"),
          .inorder = FALSE
        ) %fun% {
          lme4::glmer(
            sim.dat[, i] ~ as.factor(trt) + (1 | clust),
            family = stats::poisson(link = 'log'),
            control = lme4::glmerControl(
              optimizer = "nloptwrap",
              calc.derivs = TRUE,
              optCtrl = list(
                algorithm = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
        }
      }
      
      if (analysis == 'neg.binom') {
        if (opt == "auto") {
          mod <- lme4::glmer.nb(
            sim.dat[, 1] ~ as.factor(trt) + (1 | clust),
            control = lme4::glmerControl(
              optimizer = "nloptwrap",
              calc.derivs = TRUE,
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
          
          
        }
        my.mod <- foreach::foreach(
          i = 1:nsim,
          .options.parallel = opts,
          .packages = c("lme4", "optimx", "nloptr"),
          .inorder = FALSE
        ) %fun% {
          lme4::glmer.nb(
            sim.dat[, i] ~ as.factor(trt) + (1 | clust),
            control = lme4::glmerControl(
              optimizer = "nloptwrap",
              calc.derivs = TRUE,
              optCtrl = list(
                method = goodopt,
                starttests = FALSE,
                kkt = FALSE
              )
            )
          )
        }
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(3 / 5)
        Sys.sleep(1 / 100)
      }
      
      # re-fit models that didn't converge and
      #option to stop the function early if fits are singular
      for (i in 1:nsim) {
        converged[i] <-
          ifelse(is.null(my.mod[[i]]@optinfo$conv$lme4$messages),
                 TRUE,
                 FALSE)
      }
      singular <- which(unlist(converged) == FALSE)
      for (j in singular) {
        if (family == 'poisson') {
          my.mod[[j]] <-
            lme4::glmer(
              sim.dat[, j] ~ as.factor(trt) + (1 | clust),
              family = stats::poisson(link = 'log'),
              control = lme4::glmerControl(
                optimizer = "nloptwrap",
                calc.derivs = TRUE,
                optCtrl = list(
                  method = goodopt,
                  starttests = FALSE,
                  kkt = FALSE
                )
              )
            )
        }
        if (family == 'neg.binom') {
          my.mod[[j]] <-
            lme4::glmer.nb(
              sim.dat[, j] ~ as.factor(trt) + (1 | clust),
              control = lme4::glmerControl(
                optimizer = "nloptwrap",
                calc.derivs = TRUE,
                optCtrl = list(
                  method = goodopt,
                  starttests = FALSE,
                  kkt = FALSE
                )
              )
            )
        }
      }
      
      for (i in 1:nsim) {
        converged[i] <-
          ifelse(is.null(my.mod[[i]]@optinfo$conv$lme4$messages),
                 TRUE,
                 FALSE)
      }
      if (poor.fit.override == FALSE) {
        if (sum(isFALSE(converged), na.rm = TRUE) > (nsim * .25)) {
          stop("more than 25% of simulations are singular fit: check model specifications")
        }
      }
      
      if (!is.na(cores) & quiet == FALSE) {
        message("\r Performing null model comparisons")
      }
      # get the overall p-values (>Chisq)
      model.compare <- foreach::foreach(
        i = 1:nsim,
        .options.parallel = opts,
        .packages = "car",
        .inorder = FALSE
      ) %fun% {
        car::Anova(my.mod[[i]], type = "II")
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(4 / 5)
        Sys.sleep(1 / 100)
      }
      gc()
      # get the model summaries
      if (!is.na(cores) & quiet == FALSE) {
        message("\r Retrieving model summaries")
      }
      model.values <-
        foreach::foreach(
          i = 1:nsim,
          .options.parallel = opts,
          .packages = "car",
          .inorder = FALSE
        ) %fun% {
          summary(my.mod[[i]])
        }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(5 / 5)
        Sys.sleep(1 / 100)
      }
      # turn off parallel computing
      if (!is.na(cores)) {
        #stop the progress bar
        close(pb)
        parallel::stopCluster(cl)
      }
    } #end of GLMM method
    
    # Fit GEE (geeglm)
    if (method == 'gee') {
      sim.start <- Sys.time()
      geepack::geeglm(
        sim.dat[, 1] ~ as.factor(trt),
        family = stats::poisson(link = 'log'),
        id = clust,
        corstr = "exchangeable"
      )
      time.est = avg.iter.time * (nsim - 1) / 60
      hr.est = time.est %/% 60
      min.est = round(time.est %% 60, 3)
      
      #time limit override (for Shiny)
      if (min.est > 2 && timelimitOverride == FALSE) {
        stop(paste0(
          "Estimated completion time: ",
          hr.est,
          'Hr:',
          min.est,
          'Min'
        ))
      }
      if (quiet == FALSE) {
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
        # initialize progress bar
        if (is.na(cores)) {
          prog.bar =  progress::progress_bar$new(format = "(:spin) [:bar] :percent eta :eta",
                                                 total = 5,
                                                 clear = FALSE)
          prog.bar$tick(0)
        }
      }
      if (!is.na(cores) & quiet == FALSE) {
        message("Fitting models")
      }
      
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(2 / 5)
        Sys.sleep(1 / 100)
      }
      if (family == "poisson") {
        my.mod <- foreach::foreach(
          i = 1:nsim,
          .options.parallel = opts,
          .packages = "geepack",
          .inorder = FALSE
        ) %fun% {
          geepack::geeglm(
            sim.dat[, i] ~ as.factor(trt),
            family = stats::poisson(link = 'log'),
            id = clust,
            corstr = "exchangeable"
          )
        }
      }
      if (family == "quasipoisson") {
        my.mod <- foreach::foreach(
          i = 1:nsim,
          .options.parallel = opts,
          .packages = "geepack",
          .inorder = FALSE
        ) %fun% {
          geepack::geeglm(
            sim.dat[, i] ~ as.factor(trt),
            family = stats::quasipoisson(link = 'log'),
            id = clust,
            corstr = "exchangeable"
          )
        }
      }
      
      # check for gee convergence
      for (i in 1:length(my.mod)){
        converged[i] <- ifelse(summary(my.mod[[i]])$error == 0, TRUE, FALSE)
      }
      
      
      if (!is.na(cores) & quiet == FALSE) {
        message("Performing null model comparisons")
      }
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(3 / 5)
        Sys.sleep(1 / 100)
      }
      null.mod <- list()
      null.mod[[i]] <- geepack::geeglm(
        sim.dat[, i] ~ 1,
        family = stats::quasipoisson(link = 'log'),
        id = clust,
        corstr = "exchangeable"
      )
      # get the overall p-values (>Chisq)
      model.compare <-
        foreach::foreach(i = 1:nsim, .inorder = FALSE) %fun% {
          try(anova(my.mod[[i]], null.mod[[i]]))
        }
      if (is.na(cores) & quiet == FALSE) {
        # Iterate progress bar
        prog.bar$update(4 / 5)
        Sys.sleep(1 / 100)
      }
      if (!is.na(cores) & quiet == FALSE) {
        message("Retrieving model summaries")
      }
      # get the model summaries
      model.values <-  foreach::foreach(i = 1:nsim,
                                        .packages = "car",
                                        .inorder = FALSE) %fun% {
                                          summary(my.mod[[i]])
                                        }
      # turn off parallel computing
      if (!is.na(cores)) {
        #stop the progress bar
        close(pb)
        parallel::stopCluster(cl)
      }
    } # end of GEE method
    
    ## Output objects
    if (all.sim.data == TRUE) {
      complete.output.internal <-  list(
        "estimates" = model.values,
        "model.comparisons" = model.compare,
        "converged" = unlist(converged),
        "sim.data" = data.frame(sim.dat),
        "optimizer algorithm" = goodopt
      )
    } else {
      complete.output.internal <-  list(
        "estimates" = model.values,
        "model.comparisons" = model.compare,
        "converged" = unlist(converged),
        "optimizer algorithm" = goodopt
      )
    }
    
    return(complete.output.internal)
  } #end of function