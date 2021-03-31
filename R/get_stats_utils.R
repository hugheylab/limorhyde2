# f = function(x, coefs, nKnots, period, ...) {
#   coefs = matrix(coefs, ncol = 1)
#   b = getBasis(x, period, nKnots, intercept = TRUE)
#
#   # stopifnot('Number of basis components must match number of coefs provided' = ncol(b) == ncol(coefs))
#
#   y = b %*% coefs
#
#   return(y)
# }


# getOptimize = function(coFunc, tVec) {
#
#   y = coFunc(tVec)
#   tr = range(tVec)
#
#   oMax = optim(par = tVec[which.max(y)], fn = coFunc, lower = tr[1], upper = tr[2],
#                control = list(fnscale = -1), method = "L-BFGS-B")
#
#
#   oMin = optim(par = tVec[which.min(y)], fn = coFunc, lower = tr[1], upper = tr[2],
#                method = "L-BFGS-B")
#
#   res = data.table(peak_time = oMax$par, peak_value = oMax$value,
#                    trough_time = oMin$par, trough_value = oMin$value)
#
#
#   return(res)
#
# }

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


# getIdx = function(i, nCond, nKnots){
#
#   idx = c(i, nCond + (i-1)*nKnots + (1:nKnots))
#
#   return(idx)
#
# }

# getCond = function(mat, i, nCond, nKnots){
#
#   cIdx = getIdx(1, nCond, nKnots)
#   m0 = mat[, cIdx]
#
#   if(i == 1){
#     return(m0)
#   } else{
#
#     tIdx = getIdx(i, nCond, nKnots) #tIdx must be the same length as nKnots +1
#
#     m1 = mat[, tIdx]
#
#     m =  m0 + m1
#
#     colnames(m) = colnames(m1)
#
#     return(m) }
#
#
# }

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


# fixDiffPhase = function(x, period){
#
#   y = x %% period
#   y  = ifelse(y > period / 2, yes = y - period, no = y)
#   y  = ifelse(y < (-period / 2), yes = y + period, no = y)
#
#   return(y)
#
# }

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
