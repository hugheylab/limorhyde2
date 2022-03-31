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
#'   variable.
#' @param degree Integer indicating degree of the piecewise polynomial for the
#'   spline.
#' @param sinusoid Logical indicating whether to fit a cosinor-based model
#'   instead of a spline-based model.
#' @param timeColname String indicating the column in `metadata` containing the
#'   time at which each sample was acquired.
#' @param condColname String indicating the column in `metadata` containing the
#'   condition in which each sample was acquired. `NULL` indicates all samples
#'   came from the same condition. If not `NULL`, the model will include main
#'   effects and interactions with the terms for time.
#' @param covarColnames Character vector indicating the columns in `metadata`
#'   containing covariates to include in the model. `NULL` indicates no
#'   covariates.
#' @param sampleColname String indicating the column in `metadata` containing
#'   the name of each sample, which must correspond to the column names of `y`.
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
#' * `nKnots`: Number of knots.
#' * `degree`: As supplied above.
#' * `sinusoid`: As supplied above.
#' * `nConds`: Number of conditions.
#' * `nCovs`: Number of covariates.
#' * `lmFits`: If `keepLmFits` is `TRUE`, a list of objects from `limma` or
#'   `DESeq2`, with length equal to length of the `shifts` element.
#'
#' @eval examples1()
#'
#' @seealso [getPosteriorFit()]
#'
#' @export
getModelFit = function(
  y, metadata, period = 24, nKnots = 3L, degree = if (nKnots > 2) 3L else 2L,
  sinusoid = FALSE, timeColname = 'time', condColname = NULL,
  covarColnames = NULL, sampleColname = 'sample', nShifts = 3L,
  method = c('trend', 'voom', 'deseq2'), lmFitArgs = list(),
  eBayesArgs = if (method == 'trend') list(trend = TRUE) else list(),
  DESeqArgs = list(), keepLmFits = FALSE) {

  shift = NULL
  assertTRUE(length(dim(y)) == 2L)
  assertNames(rownames(y), type = 'unique')
  assertNames(colnames(y), type = 'unique')
  assertDataFrame(metadata)

  assertNumber(period, lower = .Machine$double.eps, finite = TRUE)

  assertFlag(sinusoid)
  if (sinusoid) {
    nKnots = 2L
    degree = 0L
  } else {
    assertNumber(nKnots, lower = 2)
    nKnots = assertCount(nKnots, coerce = TRUE)
    assertNumber(degree, lower = 1, upper = nKnots)
    degree = assertCount(degree, coerce = TRUE)}
  # after this point, nKnots should never be NULL

  assertString(sampleColname)
  assertChoice(sampleColname, colnames(metadata))
  assertNames(metadata[[sampleColname]], type = 'unique')

  assertSetEqual(colnames(y), metadata[[sampleColname]])
  y = y[, metadata[[sampleColname]]]

  assertString(timeColname)
  assertChoice(timeColname, colnames(metadata))
  assertNumeric(metadata[[timeColname]], finite = TRUE, any.missing = FALSE)

  assertString(condColname, null.ok = TRUE)
  if (!is.null(condColname)) {
    assertTRUE(condColname != timeColname)
    assertChoice(condColname, colnames(metadata))
    assert(checkCharacter(metadata[[condColname]]),
           checkFactor(metadata[[condColname]]))}

  assertCharacter(covarColnames, null.ok = TRUE)
  if (!is.null(covarColnames)) {
    assertDisjunct(covarColnames, c(timeColname, condColname))
    assertSubset(covarColnames, colnames(metadata))}

  nShifts = assertCount(nShifts, positive = TRUE, coerce = TRUE)
  method = match.arg(method)
  assertList(lmFitArgs)
  assertList(eBayesArgs)
  assertFlag(keepLmFits)

  shifts = getShifts(nShifts, nKnots, degree, period)
  m = getMetadata(metadata, timeColname, condColname, covarColnames)

  lmFits = foreach(shift = shifts) %do% {
    mShift = data.table::copy(m)
    set(mShift, j = 'time', value = mShift$time + shift)
    design = getDesign(mShift, period, nKnots, degree)

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

  sufs = rep(paste0('_shift', seq_len(length(shifts))),
             each = ncol(fit$coefficients) / length(shifts))
  colnames(fit$coefficients) = paste0(colnames(fit$coefficients), sufs)

  fit$shifts = shifts
  fit$period = period
  fit$conds = levels(m$cond) # always works

  fit$nKnots = nums[1L]
  fit$degree = degree
  fit$sinusoid = sinusoid
  fit$nConds = nums[2L]
  fit$nCovs = nums[3L]

  if (isTRUE(keepLmFits)) fit$lmFits = lmFits
  class(fit) = 'limorhyde2'
  return(fit)}
