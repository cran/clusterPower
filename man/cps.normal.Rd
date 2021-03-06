% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cps.normal.R
\name{cps.normal}
\alias{cps.normal}
\title{Power simulations for cluster-randomized trials: Parallel Designs, Normal Outcome}
\usage{
cps.normal(
  nsim = NA,
  nclusters = NA,
  nsubjects = NA,
  mu = 0,
  mu2 = NA,
  ICC = NA,
  sigma_sq = NA,
  sigma_b_sq = NA,
  ICC2 = NA,
  sigma_sq2 = NA,
  sigma_b_sq2 = NA,
  alpha = 0.05,
  method = "glmm",
  quiet = FALSE,
  allSimData = FALSE,
  seed = NA,
  poorFitOverride = FALSE,
  timelimitOverride = TRUE,
  lowPowerOverride = FALSE,
  irgtt = FALSE,
  nofit = FALSE
)
}
\arguments{
\item{nsim}{Number of datasets to simulate; accepts integer.  Required.}

\item{nclusters}{Number of clusters per condition; accepts single integer (implying equal numbers of clusters in the two groups)
or vector of length 2 (unequal number of clusters per arm).  Required.}

\item{nsubjects}{Number of subjects per cluster; accepts either a scalar (implying equal cluster sizes for the two groups), 
a vector of length two (equal cluster sizes within arm), or a vector of length \code{sum(nclusters)} 
(unequal cluster sizes within arm).  Required.}

\item{mu}{Mean in the first arm; accepts numeric, default 0.  Required..}

\item{mu2}{Mean in the second arm; accepts numeric.  Required.

At least 2 of the following must be specified:}

\item{ICC}{Intra-cluster correlation coefficient; accepts a value between 0 and 1.}

\item{sigma_sq}{Within-cluster variance; accepts numeric.}

\item{sigma_b_sq}{Between-cluster variance; accepts numeric.


The defaults for the following are all NA, implying equal variance parameters 
for the two groups. If one of the following is given, variance parameters differ 
between treatment groups, and at least 2 of the following 
must be specified:}

\item{ICC2}{Intra-cluster correlation coefficient for clusters in the second arm.}

\item{sigma_sq2}{Within-cluster variance for clusters in the second arm.}

\item{sigma_b_sq2}{Between-cluster variance for clusters in the second arm.

Optional parameters:}

\item{alpha}{Significance level; default = 0.05.}

\item{method}{Analytical method, either Generalized Linear Mixed Effects Model (GLMM, default) or 
Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee').}

\item{quiet}{When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.}

\item{allSimData}{Option to include a list of all simulated datasets in the output object.
Default = \code{FALSE}.}

\item{seed}{Option to set the seed. Default, NA, selects a seed based on the system clock.}

\item{poorFitOverride}{Option to override \code{stop()} if more than 25\% 
of fits fail to converge.}

\item{timelimitOverride}{Logical. When FALSE, stops execution if the estimated completion time
is more than 2 minutes. Defaults to TRUE.}

\item{lowPowerOverride}{Option to override \code{stop()} if the power
is less than 0.5 after the first 50 simulations and every ten simulations
thereafter. On function execution stop, the actual power is printed in the
stop message. Default = FALSE. When TRUE, this check is ignored and the
calculated power is returned regardless of value.}

\item{irgtt}{Logical. Default = FALSE. Is the experimental design an individually randomized 
group treatment trial? For details, see ?cps.irgtt.normal.}

\item{nofit}{Option to skip model fitting and analysis and instead return a dataframe with
the simulated datasets. Default = \code{FALSE}.}
}
\value{
If \code{nofit = F}, a list with the following components:
\itemize{
  \item Character string indicating total number of simulations and simulation type
  \item Number of simulations
  \item Data frame with columns "Power" (Estimated statistical power), 
               "lower.95.ci" (Lower 95\% confidence interval bound), 
               "upper.95.ci" (Upper 95\% confidence interval bound),
               "Alpha" (Probability of committing a type I or \mjseqn{\alpha} error or rejecting a true null),
               "Beta" (Probability of committing a type II error or failing to reject a false null).
               Note that non-convergent models are returned for review, 
               but not included in this calculation.
  \item Analytic method used for power estimation
  \item Significance level
  \item Vector containing user-defined cluster sizes
  \item Vector containing user-defined number of clusters in each arm
  \item Data frame reporting ICC, variance parameters, and means for each arm
  \item Vector containing expected group means based on user inputs
  \item Data frame with columns: 
                  "Estimate" (Estimate of treatment effect for a given simulation), 
                  "Std.err" (Standard error for treatment effect estimate), 
                  "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
                  "p.value", 
                  "converge", (Did the model converge?)
  \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
                  "y" (Simulated response value), 
                  "trt" (Indicator for arm), 
                  "clust" (Indicator for cluster)
                  }

If \code{nofit = T}, a data frame of the simulated data sets, containing:
\itemize{
  \item "arm" (Indicator for treatment arm)
  \item "clust" (Indicator for cluster)
  \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
  }
}
\description{
\loadmathjax

This function uses Monte Carlo methods (simulations) to estimate 
power for parallel design cluster-randomized trials with normal outcomes. Users 
can modify a variety of parameters to suit the simulations to their
desired experimental situation.

Users must specify the desired number of simulations, number of subjects per 
cluster, number of clusters per arm, expected means of the arms, and two of 
the following: ICC, within-cluster variance, or between-cluster variance.  
Defaults are provided for significance level, analytic method, progress 
updates, and whether the simulated data sets are retained.

Users have the option of specifying different variance parameters for each
arm, different numbers of clusters for each treatment group, and different numbers
of units within each cluster. 

Non-convergent models are not included in the calculation of exact confidence 
intervals.
}
\details{
The data generating model for observation \mjseqn{i} in cluster \mjseqn{j} is:
\mjsdeqn{y_{ij} \sim N(\mu + b_i, \sigma^2) }
for the first group or arm, where \mjseqn{b_i \sim N(0,\sigma_b^2)} 
, while for the second group, 
 
\mjsdeqn{y_{ij} \sim N(\mu_2 + b_i, \sigma_2^2) }
where \mjseqn{b_i \sim N(0,\sigma_{b_2}^2)}; if none of 
\mjseqn{\sigma_2^2, \sigma_{b_2}^2} or \code{ICC2} are used, then the second group uses
\mjseqn{b_i \sim N(0,\sigma_b^2)} 
and \mjseqn{y_{ij} \sim N(\mu_2 + b_i, \sigma^2)} .

All random terms are generated independent of one another.


For calls without \mjseqn{\sigma_2^2, \sigma_{b_2}^2} or \code{ICC2}, and using
\code{method="glmm"} the fitted model is:
\mjsdeqn{y_{ij}|b_i = \mu + \beta_1 x_{ij}  + b_i + e_{ij}}

with \mjseqn{\beta_1 = \mu_2 - \mu},
treatment group indicator \mjseqn{x_{ij} = 0} for the first group, with
\mjseqn{b_i \sim N(0, \sigma_b^2)} and \mjseqn{e_{ij} \sim N(0,\sigma^2)}.
In this case, both the random effects distribution and the residual distribution are the same for both
conditions.

Otherwise, for \code{method="glmm"} the fitted model is:
\mjsdeqn{y_{ij}|b_i = \mu + \beta_1 x_{ij} 
   + b_i I(x_{ij}=0) + e_{ij} I(x_{ij}=0) 
   + g_i I(x_{ij}=1) + f_{ij} I(x_{ij}=1)
   }

with \mjseqn{\beta_1}, \mjseqn{x_{ij}, b_i}, and \mjseqn{e_{ij}} as above, with
\mjseqn{g_i \sim N(0, \sigma_{b_2}^2)} and \mjseqn{f \sim N(0,\sigma_2^2)}, the 
random effects and residual distribution in the second group.
}
\section{Testing details}{
   
This function has been verified, where possible, against reference values from the NIH's GRT 
Sample Size Calculator, PASS11, \code{CRTsize::n4means}, and 
\code{clusterPower::cpa.normal}.
}

\examples{

# Estimate power for a trial with 10 clusters in each arm and 25 subjects in each 
# cluster, with an ICC of .3, sigma squared of 20 (implying sigma_b^2 of 8.57143) 
# in each group, with arm means of 1 and 4.75 in the two groups, using 100 simulated 
# data sets. The resulting estimated power should be 0.78.
   
\dontrun{

normal.sim = cps.normal(nsim = 100, nsubjects = 25, nclusters = 10, mu = 1, 
  mu2 = 4.75, ICC = 0.3, sigma_sq = 20, seed = 123)
}



# Estimate power for a trial with 5 clusters in one arm, those clusters having 25 subjects 
# each, 25 clusters in the other arm, those clusters having 5 subjects each, the first arm
# having a sigma squared of 20 and sigma_b squared of 8.57143, and the second a sigma squared
# of 9 and a sigma_b squared of 1, with estimated arm means of 1 and 4.75 in the first and 
# second groups, respectively, using 100 simulated data sets analyzed by the GEE method. 
# The estimated power should be 0.79, assuming seed = 123.

\dontrun{
normal.sim2 = cps.normal(nsim = 100, nclusters = c(5,25), nsubjects = c(25,5), mu = 1, 
  mu2 = 4.75, sigma_sq = 20,sigma_b_sq = 8.8571429, sigma_sq2 = 9, sigma_b_sq2 = 1, 
  method = "gee", seed = 123)
}


# Estimate power for a trial with 5 clusters in one arm, those clusters having
# 4, 5, 6, 7, 7, and 7 subjects each, and 10 clusters in the other arm,
# those clusters having 5 subjects each, with sigma_b_sq = .3 and and ICC of .3 in both arms.
# We have estimated arm means of 1 and 2 in the first and second arms, respectively, and we use
# 100 simulated data sets analyzed by the GLMM method.

\dontrun{
normal.sim2 = cps.normal(nsim = 100, nclusters = c(6,10), 
  nsubjects = list(c(4, 5, 6, 7, 7, 7), rep(5, times = 10)),
  mu = 1, mu2 = 2, sigma_b_sq = .3, ICC = .3, method = "glmm",
  seed = 1)
}

# The resulting estimated power (if you set seed = 1) should be about 0.76.

# Estimate power for a trial with 3 clusters in one arm, 
# those clusters having 25, 35, and 45 subjects each, and 10 clusters 
# in the other arm, those clusters having 5 subjects each, the first arm
# having a sigma squared of 20 and sigma_b squared of 8.57143, and the second a sigma squared
# of 9 and a sigma_b squared of 1, with estimated arm means of 1 and 4.75 in the first and 
# second groups, respectively, using 100 simulated data sets analyzed by the GLMM method.

\dontrun{

normal.sim2 <- cps.normal(nsim = 100, nclusters = c(3,10), 
  nsubjects = c(25, 35, 45, rep(5, times = 10)),
  mu = 1, mu2 = 4.75, sigma_sq = 20, sigma_b_sq = 8.8571429,
  sigma_sq2 = 9, sigma_b_sq2 = 1, method = "glmm")
}

# The resulting estimated power (if you set seed = 1) should be about 0.71.


}
\author{
Alexander R. Bogdan, Alexandria C. Sakrejda 
(\email{acbro0@umass.edu}), and Ken Kleinman 
(\email{ken.kleinman@gmail.com})
}
