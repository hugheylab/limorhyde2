#' @export
getPredictedValues = function(
  fit, times, coefType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  coefType = match.arg(coefType)
  if (coefType == 'posterior_mean' && is.null(fit$mashCoefficients)) {
    stop('No mash coefficients from which to calculate predicted values.')
  } else if (coefType == 'posterior_samples' && is.null(fit$mashPosteriorSamples)) {
    stop('No mash posterior samples from which to calculate predicted values.')}

  coefArray = getCoefArray(fit, coefType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  dummyColname = '_dummy' # should be distinct from other columns
  mPred = data.table(unique(times), 1)
  data.table::setnames(mPred, c(fit$timeColname, dummyColname))

  mOrigCond = unique(fit$metadata[, fit$condColname, with = FALSE])
  set(mOrigCond, j = dummyColname, value = 1)
  if (nrow(mOrigCond) > 0) {
    mPred = merge(mPred, mOrigCond, by = dummyColname, allow.cartesian = TRUE)}

  # TODO: simplify numeric covariates to overall mean
  mOrigCovar = unique(fit$metadata[, fit$covarColnames, with = FALSE])
  set(mOrigCovar, j = dummyColname, value = 1)
  if (nrow(mOrigCovar) > 0) {
    mPred = merge(mPred, mOrigCovar, by = dummyColname, allow.cartesian = TRUE)}

  mPred[, (dummyColname) := NULL]

  sampleColname = '_pred_sample' # should be distinct from other columns
  set(mPred, j = sampleColname, value = paste0('sample_', 1:nrow(mPred)))
  m = getMetadata(mPred, fit$timeColname, fit$condColname, fit$covarColnames)

  design = foreach(shift = fit$shifts, .combine = cbind) %do% {
    mNow = data.table::copy(m)
    set(mNow, j = 'time', value = mNow$time + shift)
    designNow = getDesign(mNow, fit$period, fit$nKnots)}

  predVals = foreach(postSampIdx = 1:nPostSamps, .combine = rbind) %dopar% {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)
    r = coefMat %*% t(design) / length(fit$shifts)
    colnames(r) = mPred[[sampleColname]]
    d1 = data.table::as.data.table(r, keep.rownames = 'feature')
    d2 = data.table::melt(
      d1, id.vars = 'feature', variable.name = sampleColname,
      variable.factor = FALSE)
    d3 = merge(mPred, d2, by = sampleColname, sort = FALSE)
    d3[, (sampleColname) := NULL]
    d3[, posterior_sample := postSampIdx]}

  if (nPostSamps == 1L) predVals[, posterior_sample := NULL]
  return(predVals[])}


#' @export
getPredictedIntervals = function(predVals, groupCols, mass = 0.9) {
  stopifnot(all(groupCols %in% colnames(predVals)),
            'posterior_sample' %in% colnames(predVals),
            'value' %in% colnames(predVals),
            !any(c('posterior_sample', 'value') %in% groupCols),
            length(mass) == 1L,
            is.numeric(mass),
            mass > 0,
            mass < 1)

  byCols = unique(c('feature', groupCols))
  predInts = predVals[
    , .(lower = stats::quantile(value, probs = (1 - mass) / 2),
        med = stats::quantile(value, probs = 0.5),
        upper = stats::quantile(value, probs = (1 + mass) / 2)),
    by = byCols]

  return(predInts)}
