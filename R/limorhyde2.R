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



# returns fitList
getModelFit = function(x, metadata, period, nKnots,...){

  bMat = getBasis(metadata$time, period, nKnots)
  formFull = ~ .
  design = stats::model.matrix(formFull, data = as.data.table(bMat))
  fit = limma::lmFit(x, design)
  fit = limma::eBayes(fit, trend = TRUE,...)

  return(fit)
}


getRhythmFit = function(x, metadata, timeColname = 'time', period = 24, nKnots, ...){
# things to check for
  # period is included
  stopifnot('Specify a numeric value > 0' = length(period) == 1L, is.numeric(period), period > 0)

  metadata = as.data.table(metadata)

  # set time
  setnames(metadata, timeColname, 'time')
  metadata[, time := time%%period]

  fit = getModelFit(x, metadata, period, nKnots, ...)

  return(fit)
}



getRhythmAsh = function(fit, ...){


  bHat = fit$coefficients[-1]

  # create standard error matrix
  sHat = sqrt(fit$s2.post) * fit$stdev.unscaled
  sHat = sHat[-1]

  data = mashr::mash_set_data(bHat, sHat)
  Uc = mashr::cov_canonical(data)
  resMash = mashr::mash(data,Uc,...)
  pm = get_pm(resMash)

  pm = cbind(pm, intercept = fit$coefficients[1])

  return(pm)

}



#option1
# create a function that takes as input

f = function(x, nKnots = NULL, coefs, period = 24, ...) {

  b = getBasis(x, period, nKnots)

  # stopifnot('Number of basis components must match number of coefs provided' = ncol(b) == ncol(coefs))

  y = rowSums(coefs * b)

  return(y)
}



getInitialVals = function(period = 24, nKnots = 4, coefs, step = period/1000, maximum = TRUE){

  t = seq(0, period, by = step)

  y = f(t, nKnots, coefs, period)

  if(isTRUE(maximum)){

    it = t[which.max(y)]

  } else {

    it = t[which.min(y)]

  }


  return(it)

}


getRhythmStats = function(mat, period = 24, nKnots, ...){ # just takes a mat with effects



  oMax = apply(mat, 1, function(r) optim(par = getInitialVals(nKnots = 4, coefs = r, maximum = TRUE),
                                          fn=f,
                                          control = list(fnscale = -1),
                                          method = "Brent", lower = 0, upper = period,
                                          nKnots = nKnots, coefs = r, period = period))

  # maxMat = do.call("rbind", oMax) produces list
  maxMat =  matrix(unlist(oMax), nrow= nrow(mat), byrow = TRUE)
  maxMat = maxMat[,1:2]
  colnames(maxMat) = c('xPeak', 'peak')


  opMin = apply(mat, 1, function(r) optim(par = getInitialVals(nKnots = 4, coefs = r, maximum = TRUE),
                                           fn=f,
                                           control = list(fnscale = 1),
                                           method = "Brent", lower = 0, upper = period,
                                           nKnots = nKnots, coefs = r, period = period))
  minMat =  matrix(unlist(opMin), nrow= nrow(mat), byrow = TRUE)
  minMat = minMat[,1:2]

  colnames(maxMat) = c('xTrough', 'trough')

  res = cbind(minMat, maxMat)

  return(res)

}
