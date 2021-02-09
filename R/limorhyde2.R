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


getCosinorBasis = function(x, period, intercept) {
  b = cbind(cos(x / period * 2 * pi),
            sin(x / period * 2 * pi))
  colnames(b) = c('cos', 'sin')
  b = addIntercept(b, intercept)
  return(b)}

getSplineBasis = function(x, period, nKnots, intercept) {
  knots = seq(0, period, length = nKnots + 2) #specific to b splines (must have the start, end, and middle knot values)
  b = pbs::pbs((x - min(x)) %% period, knots = knots[-c(1, length(knots))],
               Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
  colnames(b) = paste0('knot_', knots[-c(1, length(knots))])
  b = addIntercept(b, intercept)
  return(b)}


#' @export
getModelFit = function(x, md, formula, ...){

  design = stats::model.matrix(formula, data = md)
  fit = limma::lmFit(x, design)
  fit = limma::eBayes(fit, trend = TRUE,...)
  coMat = fit$coefficients
  coDT = as.data.table(coMat, keep.rownames = 'feature')
  setnames(coDT, '(Intercept)', 'intercept')

  # create standard error matrix
  seMat = sqrt(fit$s2.post) * fit$stdev.unscaled # get the coefficient standard errors.
  colnames(seMat) = paste0('se_', colnames(seMat))

  seDT = as.data.table(seMat, keep.rownames = 'feature')
  setnames(seDT, 'se_(Intercept)', 'se_intercept')


  csDT = merge(coDT, seDT, by = 'feature')

  #alternative we can just return the fit
  # fit$coefficients = coDT
  # fit$seDT = seDT
  # return(fit)

  return(csDT)
}
#' @export

getRhythmFit = function(x, metadata, timeColname = 'time', condColname = NULL, sampleColname = 'sample',
                        period = 24, nKnots = 3, ...){
# things to check for:
  # period is included
  stopifnot('Specify a numeric value > 0' = length(period) == 1L, is.numeric(period), period > 0)

  metadata = as.data.table(metadata)

  # set time and cond colnames
  setnames(metadata, timeColname, 'time')
  metadata[, time := time%%period]

  # all colnames of x found in metadata
  setnames(metadata, sampleColname, 'sample')
  stopifnot('Sample column in metadata must match expression matrix column names' =
              colnames(x) %in% metadata[, sample])

  # figure out model type
  if(is.null(nKnots)){ # if nKnots is NULL, use cosinor model

    b = getCosinorBasis(metadata$time, period, intercept = FALSE)
    colnames(b) = paste0(colnames(b), 't')
    b = data.table(b)

    formFull = ~ cost + sint
    } else{
    stopifnot(length(nKnots) == 1L, is.numeric(nKnots), nKnots > 0)
    b = getSplineBasis(metadata$time, nKnots = nKnots, period = period, intercept = FALSE)
    b = data.table(b)

    formFull = ~ .}

  #figure out design matrix
  if(is.null(condColname)){

    fitDT = getModelFit(x, b, formFull, ...)


  } else {

    # make sure cond has more than 1 condition
    setnames(metadata, condColname, 'cond')
    fitDT = foreach(condNow = unique(metadata$cond), .combine = rbind) %dopar% {
      mNow = cbind(metadata, b)
      kn = colnames(b)
      mNow = mNow[cond == condNow]
      x = x[, mNow$sample]
      fit = getModelFit(x, mNow[, ..kn], formFull, ...)
      fit[, cond := condNow]

      return(fit)}
    }
  return(coAll)
}


#' @export
getRhythmAsh = function(dt, ...){ # add mashr arguments

  #remove intercepts
  cols = colnames(dt)
  cols = cols[!(cols %like% 'intercept')]
  cols = cols[!(cols %like% 'cond')] # does not throw error if expr not found!

  # select columns for beta and s matrices
  colsCo = cols[!(cols %like% 'se')]
  colsSe = c('feature', cols[cols %like% c('se')])


  if(!('cond' %in% colnames(dt))){ # if only performing mashr once

    Bhat = as.matrix(dt[, ..colsCo], rownames = 'feature')
    Shat = as.matrix(dt[, ..colsSe], rownames = 'feature')

    data = mashr::mash_set_data(Bhat, Shat)
    U.c = mashr::cov_canonical(data)
    res.mash = mashr::mash(data,U.c,...)
    pm = get_pm(res.mash)
    colnames(pm) = paste0('pm_', colnames(pm))

    mashAll = as.data.table(pm, keep.rownames = 'feature')

    #add intercept back
    mashAll= merge(pmDT, dt[, .(feature, intercept)], by = "feature")
    # if we wanted to return raw and posterior effects, merge with original datatable
    # allDT = merge(dt, pmDT, by = "feature")
  } else { #if performing mashr multiple times ie: different conditions?

    # select columns for beta and s matrices without condition
    mashAll = foreach(condNow = unique(dt$cond), .combine = rbind) %dopar% {
      dtNow = dt[cond == condNow]
      Bhat = as.matrix(dtNow[, ..colsCo], rownames = 'feature')
      Shat = as.matrix(dtNow[, ..colsSe], rownames = 'feature')

      data = mashr::mash_set_data(Bhat, Shat)
      U.c = mashr::cov_canonical(data)
      res.mash = mashr::mash(data,U.c)
      pm = get_pm(res.mash)
      colnames(pm) = paste0('pm_', colnames(pm))
      pmDT = as.data.table(pm, keep.rownames = 'feature')

      #add intercept back
      pmDT = merge(pmDT, dt[, .(feature, intercept)], by = "feature")
      pmDT[, cond := condNow]

      return(pmDT)}
  }

  return(mashAll)

}



f = function(x, nKnots = NULL, coefs, period = 24, ...) {

  stopifnot(length(nKnots) == length(coefs))


  if(is.null(nKnots)){

    b = getCosinorBasis(x, period, intercept = FALSE)
  } else{
    stopifnot('nKnots must be a numeric value > 1' = length(nKnots) == 1L, is.numeric(nKnots), nKnots > 1)

    b = getSplineBasis(x, period, nKnots = nKnots, intercept = FALSE)}

  y = sum(coefs * b)

  return(y)
}




getInitialVals = function(period, nKnots = NULL, step = 0.01, ...){

  dt = data.table(time = seq(0, period, by = step))

  dt[, y = f(time, nKnots, period)]

  iVal = dt[max(y), x]

  return(iVal)

}


#' @export
getRhythmStats = function(dt, period, nKnots, ...){

  #






}


