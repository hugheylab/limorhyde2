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


getCoefMatOneCond = function(coefMat, condIdx, nConds, nKnots) {
  coefNow = coefMat[, c(1, (nConds + 1):(nConds + nKnots))]
  if (condIdx > 1) {
    i = nConds + (condIdx - 1) * nKnots + 1
    coefNow = coefNow + coefMat[, c(condIdx, i:(i + nKnots - 1))]}
  return(coefNow)}


getRmsAmp = function(f, co, period) {
  fSq = function(time) (f(time) - co[1L])^2
  r = sqrt(stats::integrate(fSq, 0, period)$value / period)
  return(r)}


centerCircDiff = function(x, p) {
  z = x %% p
  z = data.table::fifelse(
    z > p / 2, z - p, data.table::fifelse(z < (-p / 2), z + p, z))
  return(z)}


getCoefMatDiffCond = function(coefMat, condIdx, nConds, nKnots) {
  if (1L %in% condIdx) {
    condIdxNow = setdiff(condIdx, 1L)
    i = nConds + (condIdxNow - 1) * nKnots + 1
    coefNow = coefMat[, i:(i + nKnots - 1)]
  } else {
    i = nConds + (condIdx - 1) * nKnots + 1
    j = i + nKnots - 1
    coefNow = coefMat[, i[2]:j[2]] - coefMat[, i[1]:j[1]]}
  return(coefNow)}


getRmsDiffRhy = function(coefMat, condLevels) {
  c(period, nKnots, nConds) %<-%
    attributes(coefMat)[c('period', 'nKnots', 'nConds')]

  g = function(time) getBasis(time, period, nKnots, FALSE)
  condIdx = match(condLevels, attr(coefMat, 'condLevels'))
  coefNow = getCoefMatDiffCond(coefMat, condIdx, nConds, nKnots)

  feo = foreach(co = iterators::iter(coefNow, by = 'row'), .combine = c)

  r = feo %dopar% {
    f = function(time) (g(time) %*% t(co))^2
    rNow = sqrt(stats::integrate(f, 0, period)$value / period)}
  return(r)}
