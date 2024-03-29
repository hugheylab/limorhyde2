% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_expected.R
\name{getExpectedMeas}
\alias{getExpectedMeas}
\title{Compute expected measurements from fitted models}
\usage{
getExpectedMeas(
  fit,
  times,
  fitType = c("posterior_mean", "posterior_samples", "raw"),
  features = NULL,
  dopar = TRUE
)
}
\arguments{
\item{fit}{A 'limorhyde2' object.}

\item{times}{Numeric vector of times, in units of
\code{fit$metadata[[fit$timeColname]]}.}

\item{fitType}{String indicating which fitted models to use to compute the
expected measurements. A typical analysis using \code{limorhyde2} will be based
on 'posterior_mean', the default.}

\item{features}{Vector of names, row numbers, or logical values for
subsetting the features. \code{NULL} indicates all features.}

\item{dopar}{Logical indicating whether to run calculations in parallel if
a parallel backend is already set up, e.g., using
\code{\link[doParallel:registerDoParallel]{doParallel::registerDoParallel()}}. Recommended to minimize runtime.}
}
\value{
A \code{data.table}.
}
\description{
This function computes expected measurements (corresponding to the fitted
curves) for the specified times and features in all combinations of
conditions and covariates (if they exist).
}
\examples{
library('data.table')

y = GSE34018$y
metadata = GSE34018$metadata

fit = getModelFit(y, metadata)
fit = getPosteriorFit(fit)

measObs = mergeMeasMeta(y, metadata, features = c('13170', '12686'))
measFitMean = getExpectedMeas(
  fit, times = seq(0, 24, 0.5), features = c('13170', '12686'))
}
\seealso{
\code{\link[=getModelFit]{getModelFit()}}, \code{\link[=getPosteriorFit]{getPosteriorFit()}}, \code{\link[=getPosteriorSamples]{getPosteriorSamples()}},
\code{\link[=getExpectedMeasIntervals]{getExpectedMeasIntervals()}}
}
