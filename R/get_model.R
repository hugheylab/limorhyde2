#' Fit linear models for rhythmicity in one or more conditions
#'
#' This is the first step in an analysis using `limorhyde2`, the second is to
#' moderate the fits using [getPosteriorFit()].
#'
#' @param y Matrix-like object of measurements, with rows corresponding to
#'   features and columns to samples.
#' @param metadata data.frame containing experimental design information for
#'   each sample. Rows of `metadata` must correspond to columns of `y`. Row
#'   names are ignored.
#' @param period Number specifying the period for the time variable, in the same
#'   units as the values in the `timeColname` column.
#' @param nKnots Number of internal knots for the periodic spline for the time
#'   variable. Use `NULL` to fit a cosinor-based model instead of a spline-based
#'   model.
#' @param timeColname String indicating the column in `metadata` containing the
#'   time at which each sample was acquired.
#' @param condColname String indicating the column in `metadata` containing the
#'   condition in which each sample was acquired. `NULL` indicates all samples
#'   came from the same condition. If not `NULL`, the model will include main
#'   effects and interactions with the terms for time.
#' @param covarColnames Character vector indicating the columns in `metadata`
#'   containing covariates to include in the model. `NULL` indicates no
#'   covariates.
#' @param nShifts Number of shifted models to fit. Only used for periodic
#'   splines, not for cosinor. Do not change from the default unless you know
#'   what you're doing.
#' @param method String indicating method to estimate model coefficients. For
#'   microarray data, use 'trend'. For RNA-seq count data, use 'voom' or
#'   'deseq2'.
#' @param lmFitArgs List of arguments passed to [limma::lmFit()].
#' @param eBayesArgs List of arguments passed to [limma::eBayes()].
#' @param DESeqArgs List of arguments passed to [DESeq2::DESeq()].
#' @param keepLmFits Logical indicating whether to keep the complete fit objects
#'   from `limma` or `DESeq2`. Not needed by any functions in `limorhyde2`.
#'
#' @return A `limorhyde2` object with elements:
#'
#' * `metadata`: As supplied above, converted to a `data.table`.
#' * `timeColname`: As supplied above.
#' * `condColname`: As supplied above.
#' * `covarColnames`: As supplied above.
#' * `coefficients`: Matrix with rows corresponding to features and columns to
#'   model terms, including all shifted models.
#' * `shifts`: Numeric vector indicating amount by which timepoints were shifted
#'   for each shifted model.
#' * `period`: As supplied above.
#' * `conds`: If `condColname` is not `NULL`, a vector of unique values of
#'   the condition variable.
#' * `nKnots`: Number of knots, where 2 indicates a cosinor-based model.
#' * `nConds`: Number of conditions.
#' * `nCovs`: Number of covariates.
#' * `lmFits`: If `keepLmFits` is `TRUE`, a list of objects from `limma` or
#'   `DESeq2`, with length equal to length of the `shifts` element.
#'
#' @seealso [getPosteriorFit()]
#'
#' @export
getModelFit = function(
  y, metadata, period = 24, nKnots = 4L, timeColname = 'time',
  condColname = NULL, covarColnames = NULL, nShifts = 3L,
  method = c('trend', 'voom', 'deseq2'), lmFitArgs = list(),
  eBayesArgs = if (method == 'trend') list(trend = TRUE) else list(),
  DESeqArgs = list(), keepLmFits = FALSE) {

  shift = NULL
  assertDataFrame(metadata)
  assertTRUE(ncol(y) == nrow(metadata))

  assertNumber(period, lower = .Machine$double.eps, finite = TRUE)

  assertNumber(nKnots, lower = 3, null.ok = TRUE)
  nKnots = assertCount(nKnots, null.ok = TRUE, coerce = TRUE)
  if (is.null(nKnots)) nKnots = 2L
  # after this point, nKnots should never be NULL

  assertString(timeColname)
  assertChoice(timeColname, colnames(metadata))
  assertNumeric(metadata[[timeColname]], finite = TRUE, any.missing = FALSE)

  assertString(condColname, null.ok = TRUE)
  if (!is.null(condColname)) {
    assertTRUE(condColname != timeColname)
    assertChoice(condColname, colnames(metadata))}

  assertCharacter(covarColnames, null.ok = TRUE)
  if (!is.null(covarColnames)) {
    assertDisjunct(covarColnames, c(timeColname, condColname))
    assertSubset(covarColnames, colnames(metadata))}

  nShifts = assertCount(nShifts, positive = TRUE, coerce = TRUE)
  method = match.arg(method)
  assertList(lmFitArgs)
  assertList(eBayesArgs)
  assertFlag(keepLmFits)

  shifts = getShifts(nShifts, nKnots, period)
  m = getMetadata(metadata, timeColname, condColname, covarColnames)

  lmFits = foreach(shift = shifts) %do% {
    mShift = data.table::copy(m)
    set(mShift, j = 'time', value = mShift$time + shift)
    design = getDesign(mShift, period, nKnots)

    if (method == 'deseq2') {
      fitNow = do.call(DESeq2::DESeq, c(list(y, full = design), DESeqArgs))
    } else {
      v = if (method == 'voom') limma::voom(y, design) else y
      fitNow = do.call(limma::lmFit, c(list(v, design), lmFitArgs))
      fitNow = do.call(limma::eBayes, c(list(fitNow), eBayesArgs))}}

  fit = list(metadata = data.table::as.data.table(metadata),
             timeColname = timeColname,
             condColname = condColname,
             covarColnames = covarColnames)

  if (method == 'deseq2') {
    fit$coefficients = do.call(cbind, lapply(lmFits, stats::coef))
    fit$stdErrors = do.call(cbind, lapply(lmFits, stats::coef, SE = TRUE))
  } else {
    fit$coefficients = do.call(cbind, lapply(lmFits, `[[`, 'coefficients'))
    fit$stdErrors = do.call(
      cbind, lapply(lmFits, function(f) sqrt(f$s2.post) * f$stdev.unscaled))}

  idx = 1:(ncol(fit$coefficients) / length(shifts))
  nums = getNumKnotCondCovar(colnames(fit$coefficients)[idx])

  sufs = rep(paste0('_shift', 1:length(shifts)),
             each = ncol(fit$coefficients) / length(shifts))
  colnames(fit$coefficients) = paste0(colnames(fit$coefficients), sufs)

  fit$shifts = shifts
  fit$period = period
  fit$conds = levels(m$cond) # always works

  fit$nKnots = nums[1L]
  fit$nConds = nums[2L]
  fit$nCovs = nums[3L]

  if (isTRUE(keepLmFits)) fit$lmFits = lmFits
  class(fit) = 'limorhyde2'
  return(fit)}
