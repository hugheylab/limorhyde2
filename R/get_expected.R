#' Compute expected measurements from fitted models
#'
#' This function computes expected measurements (corresponding to the fitted
#' curves) for the specified times and features in all combinations of
#' conditions and covariates (if they exist). Register a parallel backend to
#' minimize runtime, e.g., using [doParallel::registerDoParallel()].
#'
#' @param fit A 'limorhyde2' object.
#' @param times Numeric vector of times, in units of
#'   `fit$metadata[[fit$timeColname]]`.
#' @param fitType String indicating which fitted models to use to compute the
#'   expected measurements. A typical analysis using `limorhyde2` will be based
#'   on 'posterior_mean', the default.
#' @param features Vector of names, row numbers, or logical values for
#'   subsetting the features. `NULL` indicates all features.
#'
#' @return A `data.table`.
#'
#' @seealso [getModelFit()], [getPosteriorFit()], [getPosteriorSamples()]
#'
#' @export
getExpectedMeas = function(
  fit, times, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  stopifnot(inherits(fit, 'limorhyde2'))
  fitType = match.arg(fitType)
  checkFitType(fit, fitType)

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  dummyColname = '_dummy' # should be distinct from other columns
  mNew = data.table(unique(times), 1)
  data.table::setnames(mNew, c(fit$timeColname, dummyColname))

  mOrigCond = unique(fit$metadata[, fit$condColname, with = FALSE])
  set(mOrigCond, j = dummyColname, value = 1)
  if (nrow(mOrigCond) > 0) {
    mNew = merge(mNew, mOrigCond, by = dummyColname, allow.cartesian = TRUE)}

  # TODO: simplify numeric covariates to overall mean
  mOrigCovar = unique(fit$metadata[, fit$covarColnames, with = FALSE])
  set(mOrigCovar, j = dummyColname, value = 1)
  if (nrow(mOrigCovar) > 0) {
    mNew = merge(mNew, mOrigCovar, by = dummyColname, allow.cartesian = TRUE)}

  mNew[, (dummyColname) := NULL]

  sampleColname = '_new_sample' # should be distinct from other columns
  set(mNew, j = sampleColname, value = paste0('new_sample_', 1:nrow(mNew)))
  m = getMetadata(mNew, fit$timeColname, fit$condColname, fit$covarColnames)

  design = foreach(shift = fit$shifts, .combine = cbind) %do% {
    mShift = data.table::copy(m)
    set(mShift, j = 'time', value = mShift$time + shift)
    designShift = getDesign(mShift, fit$period, fit$nKnots)}

  expectedMeas = foreach(postSampIdx = 1:nPostSamps, .combine = rbind) %dopar% {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)
    r = coefMat %*% t(design) / length(fit$shifts)
    colnames(r) = mNew[[sampleColname]]
    d1 = data.table::as.data.table(r, keep.rownames = 'feature')
    d2 = data.table::melt(
      d1, id.vars = 'feature', variable.name = sampleColname,
      variable.factor = FALSE)
    d3 = merge(mNew, d2, by = sampleColname, sort = FALSE)
    d3[, (sampleColname) := NULL]
    d3[, posterior_sample := postSampIdx]}

  if (nPostSamps == 1L) expectedMeas[, posterior_sample := NULL]
  attr(expectedMeas, 'fitType') = fitType
  return(expectedMeas[])}


#' Compute credible intervals for expected measurements
#'
#' This functions uses posterior samples to quantify uncertainty in the
#' expected measurements from fitted models.
#'
#' @param expectedMeas A `data.table` of expected measurements for posterior
#'   samples, as returned by [getExpectedMeas()].
#' @param mass Number between 0 and 1 indicating the probability mass for which
#'   to calculate the intervals.
#' @param method String indicating the type of interval: 'eti' for equal-tailed
#'   using [stats::quantile()], or 'hdi' for highest density using
#'   [HDInterval::hdi()].
#'
#' @return A `data.table` containing lower and upper bounds of the expected
#'   measurement for each combination of feature, time, and possibly condition
#'   and covariate.
#'
#' @seealso [getExpectedMeas()], [getStatsIntervals()]
#'
#' @export
getExpectedMeasIntervals = function(
  expectedMeas, mass = 0.9, method = c('eti', 'hdi')) {

  stopifnot(isTRUE(attr(expectedMeas, 'fitType') == 'posterior_samples'),
            length(mass) == 1L,
            is.numeric(mass),
            mass > 0.5,
            mass < 1)
  method = match.arg(method)

  byCols = setdiff(colnames(expectedMeas), c('posterior_sample', 'value'))
  getInterval = if (method == 'eti') getEti else getHdi
  expectedInts = expectedMeas[, getInterval(value, mass), by = byCols]

  attr(expectedInts, 'mass') = mass
  attr(expectedInts, 'method') = method
  return(expectedInts)}