#' @export
getModelFit = function(
  y, metadata, period = 24, nKnots = 4, timeColname = 'time',
  condColname = NULL, covarColnames = NULL, nShifts = 3,
  method = c('trend', 'voom'), lmFitArgs = list(),
  eBayesArgs = if (method == 'trend') list(trend = TRUE) else list()) {

  stopifnot(ncol(y) == nrow(metadata),
            is.numeric(nShifts),
            length(nShifts) == 1L)
  method = match.arg(method)

  if (is.null(nKnots) || nKnots == 2) {
    shifts = 0 # cosinor is invariant to shifts
  } else {
    knotInterval = period / (nKnots + 1)
    shiftInterval = knotInterval / nShifts
    shifts = seq(0, knotInterval - shiftInterval, shiftInterval)}

  m = getMetadata(metadata, timeColname, condColname, covarColnames)

  lmFits = foreach(shift = shifts) %do% {
    mShift = data.table::copy(m)
    set(mShift, j = 'time', value = mShift$time + shift)
    design = getDesign(mShift, period, nKnots)

    v = if (method == 'voom') limma::voom(y, design) else y
    fitNow = do.call(limma::lmFit, c(list(v, design), lmFitArgs))
    fitNow = do.call(limma::eBayes, c(list(fitNow), eBayesArgs))}

  fit = list(metadata = data.table::as.data.table(metadata),
             timeColname = timeColname,
             condColname = condColname,
             covarColnames = covarColnames,
             lmFits = lmFits)
  fit$coefficients = do.call(cbind, lapply(lmFits, `[[`, 'coefficients'))

  sufs = rep(paste0('_shift', 1:length(shifts)),
             each = ncol(fit$coefficients) / length(shifts))
  colnames(fit$coefficients) = paste0(colnames(fit$coefficients), sufs)

  fit$shifts = shifts
  fit$period = period
  fit$condLevels = levels(m$cond) # always works

  c(nK, nCon, nCov) %<-% getNumKnotCondCovar(colnames(lmFits[[1L]]))
  fit$nKnots = nK
  fit$nConds = nCon
  fit$nCovs = nCov
  return(fit)}
