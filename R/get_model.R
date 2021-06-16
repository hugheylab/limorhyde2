#' Fit a linear model to the measurements of each feature
#'
#' Given a data matrix where each row is a feature and each column is a time sample,
#' `getModelFit` returns a list object from fitting a linear model to the
#' measurements of each feature To create the linear model, the function decomposes
#' a periodic time variable into multiple linear components based on sine and
#' cosine terms or periodic spline terms of the same period.
#'
#' @param y A matrix-like data object where each row is a feature and each column
#' corresponds to a time sample.
#' @param metadata A data.table specifying experimental design information for
#' each sample. Each row is a sample with metadata given in columns.
#' @param period Number specifying the period for time variable.
#' Must be same unit as sample timepoints.
#' @param nKnots Number of knots or internal breakpoints of periodic spline.
#' @param timeColname String of column in `metadata` with the time each sample
#' was acquired.
#' @param condColname String indicating column in `metadata` with condition/group
#' name(if any) for each sample.
#' @param covarColnames String vector of covariate column names in `metadata`
#' to include in linear model.
#' @param nShifts Number of times to offset or shift time vector. Model will fit
#' data using each new shifted time vector.
#' @param method String indicating the fitting method for [limma::lmFit()].
#' Takes one of 'trend' or 'voom'.
#' @param lmFitArgs List of arguments for [limma::lmFit()].
#' @param eBayesArgs List of arguments for [limma::eBayes()].
#'
#' @return A limorhyde2 object with elements:
#'
#' * `metadata`: See `metadata` argument above.
#' * `timeColname`: See `timeColname` argument above.
#' * `condColname`: See `condColname` argument above.
#' * `covarColnames`: See `covarColnames` argument above.
#' * `lmFits`: List of linear model objects with fit results for all models.
#' * `coefficients`: A matrix with features as rows.
#' Columns are coefficient estimates for `nShift` fitted models.
#' * `shifts`: vector of shift time values
#' * `period`: See `period` argument above
#' * `condLevels`: Vector of strings indication names of conditions if available.
#' * `nKnots`: See `nKnots` argument above.
#' * `nConds`: Number of groups or conditions.
#' * `nCovs`: Number of covariates.
#'
#' @seealso [limma::lmFit()], [limma::eBayes()]
#'
#' @export
getModelFit = function(
  y, metadata, period = 24, nKnots = 4, timeColname = 'time',
  condColname = NULL, covarColnames = NULL, nShifts = 3,
  method = c('trend', 'voom'), lmFitArgs = list(),
  eBayesArgs = if (method == 'trend') list(trend = TRUE) else list(),
  keepLmFits = FALSE) {

  stopifnot(ncol(y) == nrow(metadata),
            is.numeric(nShifts),
            length(nShifts) == 1L)
  method = match.arg(method)

  shifts = getShifts(nShifts, nKnots, period)
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
             covarColnames = covarColnames)

  fit$coefficients = do.call(cbind, lapply(lmFits, `[[`, 'coefficients'))
  fit$stdErrors = do.call(
    cbind, lapply(lmFits, function(f) sqrt(f$s2.post) * f$stdev.unscaled))

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

  if (isTRUE(keepLmFits)) fit$lmFits = lmFits
  class(fit) = 'limorhyde2'
  return(fit)}
