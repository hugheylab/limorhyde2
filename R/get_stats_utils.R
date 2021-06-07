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
  c(lower, upper) %<-% range(tt)
  optMin = stats::optim(
    tt[which.min(xr)], f, method = 'L-BFGS-B', lower = lower, upper = upper)
  optMax = stats::optim(
    tt[which.max(xr)], f, method = 'L-BFGS-B', lower = lower, upper = upper,
    control = list(fnscale = -1))
  d = data.table(peak_phase = optMax$par, peak_value = optMax$value,
                 trough_phase = optMin$par, trough_value = optMin$value)
  return(d)}


getCoefMatOneCond = function(coefMat, condIdx, nConds, nKnots, nShifts) {
  nCoefs = ncol(coefMat) / nShifts
  coefKeep = foreach(j = 1:nShifts, .combine = cbind) %do% {
    coefTmp = coefMat[, (1:nCoefs) + nCoefs * (j - 1), drop = FALSE]
    coefNow = coefTmp[, c(1, (nConds + 1):(nConds + nKnots)), drop = FALSE]
    if (condIdx > 1) {
      i = nConds + (condIdx - 1) * nKnots + 1
      j = c(condIdx, i:(i + nKnots - 1))
      coefNow = coefNow + coefTmp[, j, drop = FALSE]}
    coefNow}
  return(coefKeep)}


getRmsAmp = function(f, co, period) {
  fSq = function(time) (f(time) - co[1L])^2
  r = sqrt(stats::integrate(fSq, 0, period)$value / period)
  return(r)}


centerCircDiff = function(x, p) {
  z = x %% p
  z = data.table::fifelse(
    z > p / 2, z - p, data.table::fifelse(z < (-p / 2), z + p, z))
  return(z)}


getCoefMatDiffCond = function(coefMat, condIdx, nConds, nKnots, nShifts) {
  nCoefs = ncol(coefMat) / nShifts
  coefKeep = foreach(j = 1:nShifts, .combine = cbind) %do% {
    coefTmp = coefMat[, (1:nCoefs) + nCoefs * (j - 1), drop = FALSE]
    if (1L %in% condIdx) {
      condIdxNow = setdiff(condIdx, 1L)
      i = nConds + (condIdxNow - 1) * nKnots + 1
      coefNow = coefMat[, i:(i + nKnots - 1), drop = FALSE]
    } else {
      i = nConds + (condIdx - 1) * nKnots + 1
      j = i + nKnots - 1
      coefNow = coefMat[, i[2]:j[2]] - coefMat[, i[1]:j[1], drop = FALSE]}
    coefNow}
  return(coefKeep)}


getRmsDiffRhy = function(fit, condLevels, fitType, featureIdx) {
  c(shifts, period, nKnots, nConds) %<-%
    fit[c('shifts', 'period', 'nKnots', 'nConds')]

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, FALSE)}))}

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  coefArray = coefArray[featureIdx, , , drop = FALSE]

  condIdx = match(condLevels, fit$condLevels)
  doPost = if (nPostSamps == 1L) `%do%` else `%dopar%`
  doFeat = if (nPostSamps == 1L) `%dopar%` else `%do%`

  r = doPost(foreach(postSampIdx = 1:nPostSamps, .combine = rbind), {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)
    coefNow = getCoefMatDiffCond(coefMat, condIdx, nConds, nKnots, length(shifts))

    r1 = doFeat(foreach(co = iterators::iter(coefNow, by = 'row'), .combine = c), {
      f = function(time) (g(time) %*% t(co))^2
      r2 = sqrt(stats::integrate(f, 0, period)$value / period)})

    r3 = data.table(feature = rownames(coefMat),
                    posterior_sample = postSampIdx,
                    rms_diff_rhy = r1)})

  if (nPostSamps == 1L) r[, posterior_sample := NULL]
  return(r[])}


getEti = function(v, mass) {
  r = stats::quantile(v, probs = c(1 - mass, 1 + mass) / 2)
  d = data.table(lower = r[1L], upper = r[2L])
  return(d)}


getHdi = function(v, mass) {
  r = HDInterval::hdi(v, credMass = mass)
  d = data.table(lower = r[1L], upper = r[2L])
  return(d)}
