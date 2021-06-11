#' @export
getFittedValues = function(
  fit, times, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

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
#' \code{getFittedIntervals} constructs credible intervals for the fitted values
#' of a model by group.
#'
#' @param fittedVals A data.table of posterior samples of fitted values.
#' @param groupCols A vector of column names by which posterior samples are grouped.
#' @param mass The probability mass for which to calculate the interval.
#' @param method One of 'eti' or 'hdi'. Inputting 'eti' returns an equal-tailed
#' interval (i.e., an interval with an equal probability mass in each tail);
#' inputting 'hdi' returns a highest-posterior density interval (i.e., the
#' narrowest interval containing the specified probability mass).
#'
#' @return A data.table with columns for the upper and lower bounds of the fitted
#' value for each feature in each group.
#'
#' @export
getFittedIntervals = function(
  fittedVals, groupCols, mass = 0.9, method = c('eti', 'hdi')) {

  stopifnot(isTRUE(attr(fittedVals, 'fitType') == 'posterior_samples'),
            all(groupCols %in% colnames(fittedVals)),
            !any(c('posterior_sample', 'value') %in% groupCols),
            length(mass) == 1L,
            is.numeric(mass),
            mass > 0.5,
            mass < 1)
  method = match.arg(method)

  byCols = unique(c('feature', groupCols))
  getInterval = if (method == 'eti') getEti else getHdi
  fittedInts = fittedVals[, getInterval(value, mass), by = byCols]

  attr(fittedInts, 'mass') = mass
  attr(fittedInts, 'method') = method
  return(fittedInts)}
