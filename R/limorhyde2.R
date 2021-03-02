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
  nKnots = length(cols) - nCond/nCond

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


f = function(x, coefs, nKnots, period, ...) {

  b = getBasis(x, period, nKnots)

  # stopifnot('Number of basis components must match number of coefs provided' = ncol(b) == ncol(coefs))

  y = coefs[1] + (b %*% coefs[-1])

  return(y)
}


getOptimize = function(coFunc, tVec) {

  y = coFunc(tVec)

  oMax = optim(par = tVec[which.max(y)], fn = coFunc, lower = min(tVec), upper = max(tVec),
               control = list(fnscale = -1), method = "L-BFGS-B")


  oMin = optim(par = tVec[which.min(y)], fn = coFunc, lower = min(tVec), upper = max(tVec),
               method = "L-BFGS-B")

  res = data.table(peakTime = oMax$par, peakValue = oMax$value,
                   troughTime = oMin$par, troughValue = oMin$value)


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
  t = seq(0, period, by = period/1000)

  # statsAll = foreach(codNow = (1:nCond), .combine = rbind) %do% {

  # isolate matrix for specific condition
  # maybe create function that isolates individual condition matrix by index in original matrix

  statsNow = foreach(mNow = iter(mat, by = "row"), .combine = rbind) %dopar% {

    funcR = function(x, co = mNow, p = period, nk = nKnots){

      b = getBasis(x, p, nk)

      y = co[1] + (b %*% co[-1])

      return(y)

    }


    res = getOptimize(funcR, t)
    res[, feature := rownames(mNow)]

    res[, ampl := 0.5 * (peakValue - troughValue)]

    return(res) }

  statsNow[,cond := nCond]

  return(statsNow) }