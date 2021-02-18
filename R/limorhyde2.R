#' @importFrom data.table data.table := setnames setorderv setcolorder setDT
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom zeallot %<-%
NULL


globalVariables(c('estimate', 'feature', 'variable', 'se', 'meanNow', 'j', 'ix',
                  'iy', 'jx', 'jy', 'pval_diff_amp', 'pval_diff_mesor',
                  'pval_diff_phase', 'P.Value', 'adj.P.Val', 'cond', 'varNow',
                  'qval_diff_phase', 'qval_diff_rhy', 'pval_diff_rhy',
                  '.', 'pval_rhy', 'qval_diff_amp', 'qval_diff_mesor',
                  'cond_lev', 'cond_num', 'model_comp', 'coef_idx'))




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

  fBasis = which(!(cols %like% 'cond') & cols >1) #should work with or  w/o cond
  nKnots =  length(fBasis)

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

  y = coefs[1] + rowSums(coefs[-1] * b)

  return(y)
}



getInitialVals = function(period, nKnots, coefs, step = period/1000, maximum = TRUE){


  t = seq(0, period, by = step)

  y = f(t, nKnots, coefs, period)

  if(isTRUE(maximum)){

    it = t[which.max(y)]


  } else {

    it = t[which.min(y)]

  }


  return(it)

}

getOptimize = function(coefs, period, nKnots, max) {

  initVals = getInitialVals(period, nKnots, coefs, maximum = max)

  if(isTRUE(max)) {
  o = optim(par = initVals, fn = f, lower = 0, upper = period,
            control = list(fnscale = -1), method = "L-BFGS-B",
            nKnots = nKnots, coefs = coefs, period = period)
  colN = c('xPeak','yPeak')
  } else{

    o = optim(par = initVals, fn = f, lower = 0, upper = period,
              method = "L-BFGS-B",
              nKnots = nKnots, coefs = coefs, period = period)
    colN = c('xTrough', 'yTrough')}

  res = matrix(c(o$par, o$value), nrow = 1)
  colnames(res) = colN

  return(res)


}

getRhythmStats = function(mat, period, nKnots, ...){ # just takes a mat with effects


  test = foreach(mNow = iter(mat, by = "row"), .combine = rbind) %dopar% {


    # for max
    resMax = getOptimize(coefs = mNow, period, nKnots, max = TRUE)
    #for min
    resMin = getOptimize(mNow, period, nKnots, max = FALSE)

    res = cbind(resMax, resMin)
    rownames(res) = rownames(mNow)

    res = as.data.table(res, keep.rownames = 'feature')


    return(res) }

  return(test) #return a data.table

  }
