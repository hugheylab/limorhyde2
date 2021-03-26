#' @import foreach
#' @import data.table
#' @import mashr
#' @importFrom stats optim
#' @importFrom iterators iter
#' @importFrom zeallot %<-%
#' @importFrom utils combn tail
NULL

# data.table data.table as.data.table := setnames setorderv setcolorder setDT

globalVariables(c('ampl', 'cond', 'feature', 'nKnots', 'period','mNow', 'peak_value',
                  'trough_value', 'cNum', 'condNow', 'peak_time', 'trough_time','.'))


addIntercept = function(b, intercept) {

  stopifnot('intercept argument is either TRUE or FALSE'=is.logical(intercept))

  if (intercept) {
    b = cbind(1, b)
    colnames(b)[1L] = 'intercept'}


  return(b)}


getBasis = function(x, period = 24, nKnots = 4, intercept=FALSE){

  stopifnot('Supply a numeric value > 0 for period' = length(period) == 1L,
                                                        is.numeric(period),
                                                        period > 0)

  stopifnot('Supply a numeric value > 0 for nKnots' = length(nKnots) == 1L,
            (is.numeric(nKnots) & nKnots >= 2) |is.null(nKnots))

  if(is.null(nKnots)|identical(2, nKnots)){

    b = cbind(cos(x / period * 2 * pi),
              sin(x / period * 2 * pi))
    colnames(b) = paste0('basis', 1:ncol(b))
  } else {

    knots = seq(0, period, length = nKnots + 2) # add 2 to include boundary values
    b = pbs::pbs(x %% period, knots = knots[-c(1, length(knots))],
                 Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
    colnames(b) = paste0('basis', 1:nKnots)

    b = b - (1/(nKnots+1)) }

    b = addIntercept(b, intercept)


    return(b)}

# test output w and w/o conditionsColname
getSm = function(md, timeColname, conditionsColname){

  md = as.data.table(md)
  stopifnot( 'time column is not in metadata'= length(timeColname) == 1L,
             'time column is not in metadata'= timeColname %in% colnames(md),
             is.character(timeColname))

  if(is.null(conditionsColname)){

    colsKeep = c(timeColname)
    sm = md[, colsKeep, with=FALSE]
    setnames(sm, colsKeep, 'time')
    sm[, cond := 'one']} else{

      stopifnot('Specify a string indicating the condition/treatment column name in metadata'=
                  length(conditionsColname) == 1L,
                'A column indicating the condition for each sample is not in metadata'=
                  conditionsColname %in% colnames(md),
                is.character(conditionsColname))
      colsKeep = c(timeColname, conditionsColname)

      sm = md[, colsKeep, with=FALSE] # or use :=
      setnames(sm, colsKeep, c('time','cond'))
  }

  stopifnot('time column must have numbers only' = is.numeric(sm$time))

  return(sm)
}

#' @export
getModelFit = function(x, metadata, period = 24, timeColname, conditionsColname, nKnots,...){

  sm = getSm(metadata, timeColname, conditionsColname)
  sm[, cond := as.factor(cond)]
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

  attr(fit$coefficients, 'period') = period
  attr(fit$coefficients, 'nKnots') = nKnots
  attr(fit$coefficients, 'cond') = unique(sm$cond)


  return(fit)
}

getCK = function(mat){

  cols = colnames(mat)
  nCond = which(cols == 'basis1') - 1
  nKnots = (length(cols) - nCond)/nCond

  ck = c(nCond, nKnots)

  return(ck)

}

#' @export
getRhythmAsh = function(fit, covMethod = c('canonical', 'data-driven')
  , getSigResArgs = list(), covEdArgs = list()
  , npc = attr(fit$coefficients, 'nKnots'), ...){

  c(period, nKnots, cond) %<-% attributes(fit$coefficients)[-2:-1]

  bMat = fit$coefficients
  idxRemove = getCK(bMat)[1]

  bHat = fit$coefficients[, -(1:idxRemove)]
  sHat = sqrt(fit$s2.post) * fit$stdev.unscaled[, -(1:idxRemove)]

  data = mashr::mash_set_data(bHat, sHat)
  Uc = mashr::cov_canonical(data)

  covType = match.arg(covMethod)

  if ('data-driven' %in% covMethod
    & !is.null(npc)) {

    m1by1 = mashr::mash_1by1(data)

    strong = do.call(mashr::get_significant_results
      , c(list(m1by1), getSigResArgs))

    Upca = mashr::cov_pca(data = data, subset = strong, npc = npc)

    Ued = do.call(mashr::cov_ed, c(list(data = data), list(Ulist_init = Upca)
      , list(subset = strong), covEdArgs))

  } else if ('data-driven' %in% covMethod
      & is.null(npc)) {

    stop("Data-driven method specified without argument 'npc'.")

  } else { Ued = NULL }

  resMash = mashr::mash(data,c(Uc, Ued))
  pm = resMash$result$PosteriorMean

  pm = cbind(bMat[, 1:idxRemove, drop = FALSE], pm)

  attr(pm, 'period') = period
  attr(pm, 'nKnots') = nKnots
  attr(pm, 'cond') = cond

  return(pm)}


f = function(x, coefs, nKnots, period, ...) {
  coefs = matrix(coefs, ncol = 1)
  b = getBasis(x, period, nKnots, intercept = TRUE)

  # stopifnot('Number of basis components must match number of coefs provided' = ncol(b) == ncol(coefs))

  y = b %*% coefs

  return(y)
}


getOptimize = function(coFunc, tVec) {

  y = coFunc(tVec)
  tr = range(tVec)

  oMax = optim(par = tVec[which.max(y)], fn = coFunc, lower = tr[1], upper = tr[2],
               control = list(fnscale = -1), method = "L-BFGS-B")


  oMin = optim(par = tVec[which.min(y)], fn = coFunc, lower = tr[1], upper = tr[2],
               method = "L-BFGS-B")

  res = data.table(peak_time = oMax$par, peak_value = oMax$value,
                   trough_time = oMin$par, trough_value = oMin$value)


  return(res)

}


getIdx = function(i, nCond, nKnots){

  idx = c(i, nCond + (i-1)*nKnots + (1:nKnots))

  return(idx)

}


getCond = function(mat, i, nCond, nKnots){

  cIdx = getIdx(1, nCond, nKnots)
  m0 = mat[, cIdx]

  if(i == 1){
    return(m0)
  } else{

    tIdx = getIdx(i, nCond, nKnots) #tIdx must be the same length as nKnots +1

    m1 = mat[, tIdx]

    m =  m0 + m1

    colnames(m) = colnames(m1)

    return(m) }


}


#' @export
getRhythmStats = function(mat){ # just takes a mat with effects

  c(period, nKnots, cond) %<-% attributes(mat)[-2:-1]

  t = seq(0, period, by = period/(nKnots *10))

  statsAll = foreach(cNum = 1:length(cond), .combine = rbind) %do% {

    cMat = getCond(mat, cNum, length(cond), nKnots)
    # isolate matrix for specific condition
    # maybe create function that isolates individual condition matrix by index in original matrix

    statsNow = foreach(mNow = iter(cMat, by = "row"), .combine = rbind) %dopar% {

      funcR = function(x, co = mNow, p = period, nk = nKnots){
        co2 = matrix(co, ncol = 1)
        b = getBasis(x, p, nk, intercept = TRUE)

        y = b %*% co2


        return(y)

      }


      res = getOptimize(funcR, t)
      res[, feature := rownames(mNow)]

      res[, ampl := peak_value - trough_value]

      # res = data.table(funcR(t))

      return(res) }

    if(length(cond) > 1) {

      statsNow[,cond := cond[cNum]]
    }


    return(statsNow) }

  attr(statsAll, 'period') = period
  attr(statsAll, 'nKnots') = nKnots
  attr(statsAll, 'cond') = cond

  return(statsAll)

}


fixDiffPhase = function(x, period){

  y = x %% period
  y  = ifelse(y > period / 2, yes = y - period, no = y)
  y  = ifelse(y < (-period / 2), yes = y + period, no = y)

  return(y)

}



#' @export
getDiffRhythmStats = function(dat, condIds){

  c(period, nKnots, cond) %<-% tail(attributes(dat),3)

  stopifnot(length(condIds) >=2, condIds %in% dat[, unique(cond)])

  d = dat[cond %in% condIds]
  d = d[order(factor(cond, levels = condIds))] #preserve order

  cols = setdiff(colnames(d), c('cond', 'feature'))
  d1 = d[, lapply(.SD, diff), .SDcols = cols, by  = feature]
  d1[, cond := paste(condIds, collapse = ':')]

  d1[, peak_time := fixDiffPhase(peak_time, period)]
  d1[, trough_time := fixDiffPhase(trough_time, period)]

  setnames(d1, cols, paste0("diff_", cols))


  return(d1)

}



