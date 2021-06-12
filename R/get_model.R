#' fit a linear model to the measurements of each gene
#'
#' Given a data matrix where each row is a gene and each column is a time sample,
#'  `getModelFit` returns a list object from fitting a linear model to the
#'  expression of each gene. To create the linear model, the function decomposes
#'  a periodic time variable into multiple linear components based on sine and
#'  cosine terms or periodic spline terms of the same period.
#'
#' @param y a matrix-like data object where each row is a gene and each column
#' corresponds to a time sample.
#' @param metadata a data.table specifying experimental design information for
#' each sample. Each row is a sample and metadata given in columns.
#' @param period number specifying the period for periodic time variable.
#' Must be same unit as sample timepoints.
#' @param nKnots number of knots or internal breakpoints of periodic spline
#' @param timeColname string of column in `metadata` with times samples
#' were acquired.
#' @param condColname string of column in `metadata` with condition/group
#' names(if any) for each sample
#' @param covarColnames string vector of covariate column names in `metadata`
#' to include in linear model
#' @param nShifts number of times to offset or shift time vector. Model will fit
#' data using each new shifted time vector.
#' @param method string indicating the fitting method for limma `mFit`.
#' Takes one of `‘trend’ and ‘voom’`
#' @param lmFitArgs list of arguments for limma `lmFit`
#' @param eBayesArgs list of arguments for `limma::eBayes`
#'
#' @return a `LimoRhyde2` class object with the results of
#' `limma::lmFit` including:
#'
#' * `coefficients` a matrix with rows for each feature.
#' Columns are coefficient estimates for `nShift` fitted models.
#'
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
  class(fit) = 'limorhyde2'
  return(fit)}
