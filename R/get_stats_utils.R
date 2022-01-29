checkFitType = function(fit, fitType) {
  if (fitType == 'posterior_mean' && is.null(fit$mashCoefficients)) {
    stop('No posterior mean, please run getPosteriorFit.')
  } else if (fitType == 'posterior_samples' && is.null(fit$mashPosteriorSamples)) {
    stop('No posterior samples, please run getPosteriorSamples.')}
  invisible()}


getCoefArray = function(fit, fitType) {
  coefName = switch(fitType,
                    raw = 'coefficients',
                    posterior_mean = 'mashCoefficients',
                    posterior_samples = 'mashPosteriorSamples')
  coefArray = fit[[coefName]]

  if (fitType != 'posterior_samples') {
    coefArray = array(coefArray, dim = c(dim(coefArray), 1L),
                      dimnames = dimnames(coefArray))}
  return(coefArray)}


getOptima = function(f, tt) {
  xr = f(tt)
  r = range(tt)
  optMin = stats::optim(
    tt[which.min(xr)], f, method = 'L-BFGS-B', lower = r[1L], upper = r[2L])
  optMax = stats::optim(
    tt[which.max(xr)], f, method = 'L-BFGS-B', lower = r[1L], upper = r[2L],
    control = list(fnscale = -1))
  d = data.table(peak_phase = optMax$par, peak_value = optMax$value,
                 trough_phase = optMin$par, trough_value = optMin$value)
  return(d)}


getShiftCols = function(shiftIdx, nCoefs) {
  x = (1:nCoefs) + nCoefs * (shiftIdx - 1)
  return(x)}


getBasisCols = function(condIdx, nConds, nKnots) {
  x = nConds + (1 + (condIdx - 1) * nKnots):(condIdx * nKnots)
  return(x)}


getCoefMatOneCond = function(coefMat, condIdx, nConds, nKnots, nShifts) {
  j = NULL
  coefKeep = foreach(j = 1:nShifts, .combine = cbind) %do% {
    shiftCols = getShiftCols(j, ncol(coefMat) / nShifts)
    coefTmp = coefMat[, shiftCols, drop = FALSE]
    i = getBasisCols(condIdx, nConds, nKnots)
    coefNow = coefTmp[, c(condIdx, i), drop = FALSE]
    if (condIdx > 1) coefNow[, 1L] = coefNow[, 1L] + coefTmp[, 1L]
    coefNow}
  return(coefKeep)}


getRmsAmp = function(f, co, period) {
  fSq = function(time) (f(time) - co[1L])^2
  r = sqrt(stats::integrate(fSq, 0, period)$value / period)
  return(r)}


circMean = function(amp, phase, period) {
  w = 2 * pi * phase / period
  x = amp * cos(w)
  y = amp * sin(w)
  mx = mean(x)
  my = mean(y)
  mamp = sqrt(mx^2 + my^2)
  mphase = atan2(my, mx) * period / (2 * pi)
  return(c(mamp, mphase))}


centerCircDiff = function(x, p) {
  z = x %% p
  z = data.table::fifelse(
    z > p / 2, z - p, data.table::fifelse(z < (-p / 2), z + p, z))
  return(z)}


getCoefMatDiffCond = function(coefMat, condIdx, nConds, nKnots, nShifts) {
  j = k = NULL
  coefKeep = foreach(j = 1:nShifts, .combine = cbind) %do% {
    shiftCols = getShiftCols(j, ncol(coefMat) / nShifts)
    coefTmp = coefMat[, shiftCols, drop = FALSE]
    i = foreach(k = 1:2, .combine = rbind) %do% {
      getBasisCols(condIdx[k], nConds, nKnots)}
    coefNow = coefTmp[, i[2L, ], drop = FALSE] - coefTmp[, i[1L, ], drop = FALSE]
    coefNow}
  return(coefKeep)}


getRmsDiffRhy = function(fit, conds, fitType, featureIdx, dopar) {
  shifts = period = nKnots = nConds = period = postSampIdx = co = NULL
  c(shifts, period, nKnots, nConds) %<-%
    fit[c('shifts', 'period', 'nKnots', 'nConds')]

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, FALSE)}))}

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  coefArray = coefArray[featureIdx, , , drop = FALSE]

  condIdx = match(conds, fit$conds)
  doPost = if (nPostSamps == 1L | !dopar) `%do%` else `%dopar%`
  doFeat = if (nPostSamps == 1L & dopar) `%dopar%` else `%do%`

  r = doPost(foreach(postSampIdx = 1:nPostSamps, .combine = rbind), {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)
    coefNow = getCoefMatDiffCond(coefMat, condIdx, nConds, nKnots, length(shifts))

    r1 = doFeat(foreach(co = iterators::iter(coefNow, by = 'row'), .combine = c), {
      f = function(time) (g(time) %*% t(co))^2
      r2 = sqrt(stats::integrate(f, 0, period)$value / period)})

    r3 = data.table(feature = rownames(coefMat),
                    posterior_sample = postSampIdx,
                    rms_diff_rhy = r1)})

  if (nPostSamps == 1L) set(r, j = 'posterior_sample', value = NULL)
  return(r)}


getEti = function(v, mass) {
  r = stats::quantile(v, probs = c(1 - mass, 1 + mass) / 2)
  d = data.table(lower = r[1L], upper = r[2L])
  return(d)}


getHdi = function(v, mass) {
  r = HDInterval::hdi(v, credMass = mass)
  d = data.table(lower = r[1L], upper = r[2L])
  return(d)}


getSignedAmp = function(amp, phase, period) {
  meanPhase = circMean(amp, phase, period)[2L]
  idxFlip = abs(centerCircDiff(phase - meanPhase, period)) > period / 4
  ampSigned = amp * (1 - 2 * idxFlip)
  return(ampSigned)}
