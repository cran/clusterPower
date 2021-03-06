% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpa.sw.binary.R
\name{cpa.sw.binary}
\alias{cpa.sw.binary}
\title{Power simulations for cluster-randomized trials: Stepped Wedge Design, Binary Outcome}
\usage{
cpa.sw.binary(
  nclusters = NA,
  steps = NA,
  nsubjects = NA,
  alpha = 0.05,
  timeEffect = 0,
  ICC = NA,
  p0 = NA,
  p1 = NA,
  tol = 1e-05,
  GQ = 100,
  quiet = FALSE
)
}
\arguments{
\item{nclusters}{Number of clusters; accepts non-negative integer scalar (required).}

\item{steps}{Number of crossover steps; Accepts positive scalar indicating the total
number of steps, NOT including the baseline (required).}

\item{nsubjects}{Number of subjects per cluster; accepts a scalar. Equal cluster sizes
are assumed (required).}

\item{alpha}{Significance level (default=0.05).}

\item{timeEffect}{Expected time effect over the entire study period (assumed to be linear
across time steps); accepts numeric (required). Default = 0 (no time effects).}

\item{ICC}{Intracluster correlation coefficient as defined by Hussey and Hughes (2007) 
for participants at first time step; accepts numeric (required).}

\item{p0}{Estimated baseline effect; accepts numeric (required).}

\item{p1}{Estimated treatment effect; accepts numeric (required).}

\item{tol}{Machine tolerance. Accepts numeric. Default is 1e-5.}

\item{GQ}{Number of quadriture points used in Gaussian Legendre integration; accepts
a scalar. Default is 100.}

\item{quiet}{Suppresses the progress bar; logical. Default is FALSE.}
}
\value{
The estimated power.
}
\description{
This function uses a modified Cox method to determine power for stepped
wedge cluster-randomized controlled trials. Users can modify a variety of
parameters to suit their desired experimental situation.
}
\details{
The stepped wedge trial design is a type of cross-over
design in which clusters change treatments in waves. Initially all the
clusters receive the same standard treatment, and at the end of the trial all
of the clusters will be receiving the treatment of interest. More than one
cluster can change treatments in a wave, but the order in which clusters
change treatments is randomly determined. The outcome of interest is assessed
in each cluster during each wave.

Users must specify the number of subjects per cluster, number of clusters,
the number of time steps, the baseline effect, the expected treatment effect,
expected absolute difference between treatment
arms, ICC, and time effect.
}
\note{
Much of the FORTRAN code for this package was kindly provided by Dr. Zhou.
}
\examples{

# Estimate power for a trial with 3 steps and 9 clusters at the 
# initiation of the study. Those 
# clusters have 14 subjects each with no time effects. 
# We estimated arm outcome proportions of 
# 0.2 (pre-treatment) and 0.31 (post-treatment) and intracluster 
# correlation coefficient (ICC) of 0.05. 
# The resulting power should be 0.7992842.

\dontrun{
sw.bin <- cpa.sw.binary(nclusters = 9,
  steps = 3,
  nsubjects = 14,
  timeEffect = 0,
  ICC = 0.05,
  p1 = 0.31,
  p0 = 0.2)
}

}
\references{
Zhou X, Liao X, Kunz LM, Normand ST, Wang M, Spiegelman D. A maximum
likelihood approach to power calculations for stepped wedge designs of binary
outcomes. Biostatistics. 2020 Jan 1;21(1):102-121. doi: 10.1093/biostatistics/kxy031

Hussey, MA AND Hughes, JP. (2007). Design and analysis of stepped wedge 
cluster randomized trials. Contemporary Clinical Trials 28, 182–191.
}
\author{
Alexandria C. Sakrejda (\email{acbro0@umass.edu})

Ken Kleinman (\email{ken.kleinman@gmail.com})
}
