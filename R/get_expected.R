#' Compute expected measurements from fitted models
#'
#' This function computes expected measurements (corresponding to the fitted
#' curves) for the specified times and features in all combinations of
#' conditions and covariates (if they exist).
#'
#' @param fit A 'limorhyde2' object.
#' @param times Numeric vector of times, in units of
#'   `fit$metadata[[fit$timeColname]]`.
#' @param fitType String indicating which fitted models to use to compute the
#'   expected measurements. A typical analysis using `limorhyde2` will be based
#'   on 'posterior_mean', the default.
#' @param features Vector of names, row numbers, or logical values for
#'   subsetting the features. `NULL` indicates all features.
#' @param dopar Logical indicating whether to run calculations in parallel if
#'   a parallel backend is already set up, e.g., using
#'   [doParallel::registerDoParallel()]. Recommended to minimize runtime.
#'
#' @return A `data.table`.
#'
#' @eval examples3()
#'
#' @seealso [getModelFit()], [getPosteriorFit()], [getPosteriorSamples()],
#'   [getExpectedMeasIntervals()]
#'
#' @export
getExpectedMeas = function(
  fit, times, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL, dopar = TRUE) {

  shift = postSampIdx = posterior_sample = NULL

  assertClass(fit, 'limorhyde2')
  assertNumeric(times, finite = TRUE, any.missing = FALSE)
  fitType = match.arg(fitType)
  checkFitType(fit, fitType)
  assertFlag(dopar)

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  dummyColname = '_dummy' # should be distinct from other columns
  mNew = data.table(unique(times), 1)
  data.table::setnames(mNew, c(fit$timeColname, dummyColname))

  mOrigCond = unique(fit$metadata[, fit$condColname, with = FALSE])
  if (nrow(mOrigCond) > 0) {
    set(mOrigCond, j = dummyColname, value = 1)
    mNew = merge(mNew, mOrigCond, by = dummyColname, allow.cartesian = TRUE)}

  # TODO: simplify numeric covariates to overall mean
  mOrigCovar = unique(fit$metadata[, fit$covarColnames, with = FALSE])
  if (nrow(mOrigCovar) > 0) {
    set(mOrigCovar, j = dummyColname, value = 1)
    mNew = merge(mNew, mOrigCovar, by = dummyColname, allow.cartesian = TRUE)}

  mNew[, (dummyColname) := NULL]

  sampleColname = '_new_sample' # should be distinct from other columns
  set(mNew, j = sampleColname, value = paste0('new_sample_', seq_len(nrow(mNew))))
  m = getMetadata(mNew, fit$timeColname, fit$condColname, fit$covarColnames)

  design = foreach(shift = fit$shifts, .combine = cbind) %do% {
    mShift = data.table::copy(m)
    set(mShift, j = 'time', value = mShift$time + shift)
    designShift = getDesign(mShift, fit$period, fit$nKnots, fit$degree)}

  reg = foreach::getDoParRegistered()
  doOp = if (dopar && reg) `%dopar%` else `%do%`

  expectedMeas = doOp(foreach(postSampIdx = 1:nPostSamps, .combine = rbind), {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)
    r = coefMat %*% t(design) / length(fit$shifts)
    colnames(r) = mNew[[sampleColname]]
    d1 = data.table::as.data.table(r, keep.rownames = 'feature')
    d2 = data.table::melt(
      d1, id.vars = 'feature', variable.name = sampleColname,
      variable.factor = FALSE)
    d3 = merge(mNew, d2, by = sampleColname, sort = FALSE)
    set(d3, j = sampleColname, value = NULL)
    set(d3, j = 'posterior_sample', value = postSampIdx)})

  if (nPostSamps == 1L) expectedMeas[, posterior_sample := NULL]
  setattr(expectedMeas, 'fitType', fitType)
  return(expectedMeas)}


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
#' @eval examples4()
#'
#' @seealso [getExpectedMeas()], [getStatsIntervals()]
#'
#' @export
getExpectedMeasIntervals = function(
  expectedMeas, mass = 0.9, method = c('eti', 'hdi')) {

  value = NULL
  assertDataTable(expectedMeas)
  assertTRUE(attr(expectedMeas, 'fitType') == 'posterior_samples')
  assertNumber(mass, lower = 0.5, upper = 1 - .Machine$double.eps)
  method = match.arg(method)

  byCols = setdiff(colnames(expectedMeas), c('posterior_sample', 'value'))
  getInterval = if (method == 'eti') getEti else getHdi
  expectedInts = expectedMeas[, getInterval(value, mass), by = byCols]

  setattr(expectedInts, 'mass', mass)
  setattr(expectedInts, 'method', method)
  return(expectedInts)}


#' Merge measurements and metadata
#'
#' This function is useful for plotting time-courses for individual features.
#'
#' @param y Matrix-like object of measurements, with rows corresponding to
#'   features and columns to samples.
#' @param metadata data.frame containing experimental design information for
#'   each sample. Rows of `metadata` must correspond to columns of `y`. Row
#'   names are ignored.
#' @param features Vector of names, row numbers, or logical values for
#'   subsetting the features. `NULL` indicates all features.
#' @param sampleColname String indicating the column in `metadata` containing
#'   the name of each sample, which must correspond to the column names of `y`.
#'
#' @return A `data.table` with one row for each sample-feature pair.
#'
#' @eval examples3()
#'
#' @seealso [getExpectedMeas()]
#'
#' @export
mergeMeasMeta = function(y, metadata, features = NULL, sampleColname = 'sample') {
  assertDataFrame(metadata)
  assertTRUE(ncol(y) == nrow(metadata))
  assertDisjunct(sampleColname, c('feature', 'meas'))
  assertChoice(sampleColname, colnames(metadata))

  if (!is.matrix(y)) y = as.matrix(y) # in case y is a DGEList or DESeq object
  if (!is.null(features)) y = y[features, , drop = FALSE]

  d = data.table(y, keep.rownames = 'feature')
  d = data.table::melt(
    d, id.vars = 'feature', variable.name = sampleColname, value.name = 'meas',
    variable.factor = FALSE)
  d = merge(metadata, d, by = sampleColname)
  return(d)}
