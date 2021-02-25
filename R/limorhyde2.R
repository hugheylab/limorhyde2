#' @importFrom data.table data.table := setnames setorderv setcolorder setDT
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom zeallot %<-%
NULL




addIntercept = function(b, intercept) {
  if (intercept) {
    b = cbind(1, b)
    colnames(b)[1L] = 'intercept'}
  return(b)}


getBasis = function(x, period = 24, nKnots = 4,intercept=FALSE){

  if(is.null(nKnots)){

    b = cbind(cos(x / period * 2 * pi),
              sin(x / period * 2 * pi))
    colnames(b) = paste0('basis', 1:ncol(b))
  } else {

    knots = seq(0, period, length = nKnots + 2) # add 2 to include boundary values
    b = pbs::pbs(x %% period, knots = knots[-c(1, length(knots))],
                 Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
    colnames(b) = paste0('basis', 1:nKnots)}

    b = addIntercept(b, intercept)

    return(b)}

getSm = function(md, timeColname, conditionsColname){

  md = as.data.table(md)
  colsKeep = c(timeColname, conditionsColname)
  sm = md[, ..colsKeep]

  if(is.null(conditionsColname)){
  setnames(sm, colsKeep, c('time'))
  } else{

    setnames(sm, colsKeep, c('time', 'cond'))
  }

  return(sm)
}

# returns fitList
getModelFit = function(x, metadata, period = 24, timeColname, conditionsColname, nKnots,...){

  sm = getSm(metadata, timeColname, conditionsColname)
  # things to check for
  # period is included
  stopifnot('Specify a numeric value > 0' = length(period) == 1L, is.numeric(period), period > 0)



  bMat = getBasis(sm$time, period, nKnots)

  if(is.null(conditionsColname)){
    bMat = as.data.table(bMat)
    formFull = ~ .

    } else{
      bMat = cbind(sm[, .(cond)], bMat)
      formFull = ~ cond * .}

  design = stats::model.matrix(formFull, data = bMat)

  fit = limma::lmFit(x, design)
  fit = limma::eBayes(fit, trend = TRUE,...)


  return(fit)
}


getCK = function(mat){

  cols = colnames(mat)
  nCond = which(cols == 'basis1') - 1
  nKnots =  (length(colnames) - nCond)/nCond

  ck = c(nCond, nKnots)

  return(ck)

}


getRhythmAsh = function(fit, ...){

  bMat = fit$coefficients
  idxRemove = getCK(bMat)[1]

  bHat = fit$coefficients[, -(1:idxRemove)]
  sHat = sqrt(fit$s2.post) * fit$stdev.unscaled[, -(1:idxRemove)]

  data = mashr::mash_set_data(bHat, sHat)
  Uc = mashr::cov_canonical(data)
  resMash = mashr::mash(data,Uc)
  pm = get_pm(resMash)

  pm = cbind(bMat[, 1:idxRemove, drop = FALSE], pm)

  return(pm)

}

f = function(x, nKnots, coefs, period, ...) {

  b = getBasis(x, period, nKnots)

  # stopifnot('Number of basis components must match number of coefs provided' = ncol(b) == ncol(coefs))

  y = coefs[1] + (b %*% coefs[-1])

  return(y)
}



getInitialVals = function(coefs, period, nKnots, step = period/1000){

  t = seq(0, period, by = step)

  y = f(t, nKnots, coefs, period)

  tInit = c(t[which.max(y)], t[which.min(y)], period)

  return(tInit)

}


getOptimize = function(initVals) {

  period = initVals[3]
  oMax = optim(par = initVals[1], fn = f, lower = 0, upper = period,
               control = list(fnscale = -1), method = "L-BFGS-B")

  oMin = optim(par = initVals[2], fn = f, lower = 0, upper = period, method = "L-BFGS-B")

  res = data.table(xPeak = oMax$par, yPeak = oMax$value,
                   xTrough = oMin$par, yTrough = oMin$value)

  return(res)

}
getDiffRhythmStats = function(mat){

  ck = getCK(mat)
  nCond = ck[1]
  nKnots = ck[2]
  cntrlIdx = c(1, nCond + 1:nKnots)
  controlMat = mat[, cntrlIdx]

  for(i in 2:nCond){

    idx = c(i, (1:nKnots)+i*nKnots)
    matNow = mat[, idx]
    mat[, idx] = controlMat + matNow

  }
  return(mat)

}
getRhythmStats = function(mat, period){ # just takes a mat with effects

  ck = getCK(mat)
  nCond = ck[1]
  nKnots = ck[2]

  statsAll = foreach(codNow = 1:nCond, .combine = rbind) %dopar% {

    matCond = getDiffRhythmStats(mat)

    statsNow = foreach(mNow = iter(mat, by = "row"), .combine = rbind) %dopar% {

      vals = getInitialVals(mNow, period, nKnots)
      res = getOptimize(vals)
      res[, feature := rownames(mNow)]

      return(res) }

    return(statsNow) }

  return(statsAll)

}
