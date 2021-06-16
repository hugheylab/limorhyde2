#' Compute the predicted values from a fitted model
#'
#' Calculates the fitted values for selected features at given times for all combinations of conditions and covariates.
#'
#' @param fit A 'limorhyde2' object, as provided by `getModelFit` or `getPosteriorFit`.
#' @param times A vector of timepoints at which to compute fitted response values
#' @param fitType String indicating whether to calculate statistics on the posterior mean, posterior samples, or raw model fit.
#' @param features a vector of selected feature names, row numbers, or logical values for which to calculate fitted values
#'
#' @return a data.table with fitted values for each feature (and possibly feature-condition or feature-covariate pairs) at each of the specified timepoints.
#'
#' @export
getFittedValues = function(
  fit, times, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  stopifnot(inherits(fit, 'limorhyde2'))

  fitType = match.arg(fitType)
  if (fitType == 'posterior_mean' && is.null(fit$mashCoefficients)) {
    stop('No posterior mean to calculate fitted values, please run getPosteriorFit.')
  } else if (fitType == 'posterior_samples' && is.null(fit$mashPosteriorSamples)) {
    stop('No posterior samples to calculate fitted values, please run getPosteriorSamples.')}

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

  fittedVals = foreach(postSampIdx = 1:nPostSamps, .combine = rbind) %dopar% {
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

  if (nPostSamps == 1L) fittedVals[, posterior_sample := NULL]
  attr(fittedVals, 'fitType') = fitType
  return(fittedVals[])}


#' Get credible intervals for fitted values
#'
#' `getFittedIntervals` constructs credible intervals for the fitted values of
#' a model by group.
#'
#' @param fittedVals A data.table of posterior samples for fitted values from `getFittedValues`.
#' @param mass The probability mass for which to calculate the interval.
#' @param method String for type of interval: 'eti' for equal-tailed or 'hdi' for highest (posterior) density. Equal-tailed intervals use [stats::quantiles], while HPD intervals use [HDInterval::hdi].
#'
#' @return A data.table with columns for the upper and lower bounds of the fitted value for each feature (and possibly feature-condition or feature-covariate combination).
#'
#' @seealso [getFittedValues], [stats::quantiles], [HDInterval::hdi]
#'
#' @export
getFittedIntervals = function(
  fittedVals, mass = 0.9, method = c('eti', 'hdi')) {

  stopifnot(isTRUE(attr(fittedVals, 'fitType') == 'posterior_samples'),
            length(mass) == 1L,
            is.numeric(mass),
            mass > 0.5,
            mass < 1)
  method = match.arg(method)

  byCols = setdiff(colnames(fittedVals), c('posterior_sample', 'value'))
  getInterval = if (method == 'eti') getEti else getHdi
  fittedInts = fittedVals[, getInterval(value, mass), by = byCols]

  attr(fittedInts, 'mass') = mass
  attr(fittedInts, 'method') = method
  return(fittedInts)}
