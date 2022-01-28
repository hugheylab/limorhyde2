getMash = function(coefs, ses, covMethod, getSigResArgs, npc, covEdArgs, ...) {
  md = mashr::mash_set_data(coefs, ses)
  uc = if (covMethod == 'data-driven') NULL else mashr::cov_canonical(md)

  if (covMethod == 'canonical') {
    ued = NULL
  } else {
    m1by1 = mashr::mash_1by1(md)
    strong = do.call(mashr::get_significant_results, c(list(m1by1), getSigResArgs))
    upca = mashr::cov_pca(md, npc = min(npc, length(strong)), subset = strong)
    ued = do.call(mashr::cov_ed, c(list(md, upca, strong), covEdArgs))}

  mc = mashr::mash(md, c(uc, ued), ...)
  return(list(mashData = md, mashFit = mc))}


#' Compute posterior fit for linear models for rhythmicity
#'
#' This is the second step in an analysis using `limorhyde2`, the first is to
#' fit linear models using [getModelFit()]. This function obtains posterior
#' estimates of coefficients using multivariate adaptive shrinkage (mash), which
#' learns patterns in the data and accounts for noise in the original fits. The
#' defaults for arguments should work well in most cases, so only change them if
#' you know what you're doing.
#'
#' @param fit A `limorhyde2` object.
#' @param covMethod String indicating the type(s) of covariance matrices to use
#'   for the mash fit.
#' @param getSigResArgs List of arguments passed to
#'   [mashr::get_significant_results()]. Only used if `covMethod` is
#'   'data-driven' or 'both'.
#' @param npc Number of principal components passed to [mashr::cov_pca()]. Only
#'   used if `covMethod` is 'data-driven' or 'both'.
#' @param covEdArgs List of arguments passed to [mashr::cov_ed()]. Only used if
#'   `covMethod` is 'data-driven' or 'both'.
#' @param overwrite Logical for whether to recompute the mash fit if it already
#'   exists.
#' @param ... Additional arguments passed to [mashr::mash()].
#'
#' @return A `limorhyde2` object containing everything in `fit` with added or
#'   updated elements:
#'
#' * `mashData`: list of `mash` data objects
#' * `mashFits`: list of `mash` fit objects
#' * `mashCoefficients`: Matrix of posterior mean coefficients, with rows
#'   corresponding to features and columns to model terms.
#' * `mashIdx`: Vector indicating which model terms were included in the mash
#'   fit.
#'
#' @seealso [getModelFit()], [getRhythmStats()], [getExpectedMeas()]
#'
#' @export
getPosteriorFit = function(
  fit, covMethod = c('data-driven', 'canonical', 'both'), getSigResArgs = list(),
  npc = fit$nKnots, covEdArgs = list(), overwrite = FALSE, ...) {

  shifts = nKnots = nConds = NULL

  assertClass(fit, 'limorhyde2')
  covMethod = match.arg(covMethod)
  assertList(getSigResArgs)
  assertCount(npc, positive = TRUE)
  assertList(covEdArgs)
  assertFlag(overwrite)
  assertTRUE(overwrite || is.null(fit$mashFits))

  co = fit$coefficients
  c(shifts, nKnots, nConds) %<-% fit[c('shifts', 'nKnots', 'nConds')]

  nShifts = length(shifts)
  nCoefs = ncol(co) / nShifts
  mashData = list()
  mashFits = list()
  mashIdx = c()
  i = 1L

  # TODO: what if nConds == 2, so only one condition coef, just run ash?
  # std errors for mean_value are so small, shrinkage has very little effect
  # if (nConds > 1) {
  #   idx = rep(2:nConds, nShifts) +
  #     rep(seq(0, nShifts * nCoefs - 1, nCoefs), each = nConds - 1)
  #   mashes[[i]] = getMash(co[, idx], fit$stdErrors[, idx], covMethod,
  #                          getSigResArgs, npc, covEdArgs, ...)
  #   mashIdx = c(mashIdx, idx)
  #   i = i + 1L}

  for (condIdx in 1:nConds) {
    idx = rep(getBasisCols(condIdx, nConds, nKnots), nShifts) +
      rep(seq(0, nShifts * nCoefs - 1, nCoefs), each = nKnots)
    mashOutput = getMash(co[, idx], fit$stdErrors[, idx], covMethod,
                         getSigResArgs, npc, covEdArgs, ...)
    mashData[[i]] = mashOutput$mashData
    mashFits[[i]] = mashOutput$mashFit
    mashIdx = c(mashIdx, idx)
    i = i + 1L}

  co[, mashIdx] = do.call(cbind, lapply(mashFits, ashr::get_pm))

  fit$mashData = mashData
  fit$mashFits = mashFits
  fit$mashCoefficients = co
  fit$mashIdx = mashIdx
  return(fit)}


#' Draw samples from posterior distributions of fitted models
#'
#' This is an optional step in an analysis using `limorhyde2`, and is useful for
#' quantifying uncertainty in posterior estimates of fitted curves and rhythmic
#' statistics. The function calls [mashr::mash_compute_posterior_matrices()].
#'
#' @param fit A `limorhyde2' object containing posterior fits.
#' @param nPosteriorSamples Number of samples to draw from each posterior
#'   distribution.
#' @param overwrite Logical indicating whether to recompute posterior samples if
#'   they already exist.
#'
#' @return A `limorhyde2` object containing everything in `fit` with added or
#'   updated element:
#'
#' * `mashPosteriorSamples`: a three-dimensional array of coefficients, with dim
#'   1 corresponding to features, dim 2 to model terms, and dim 3 to posterior
#'   samples.
#'
#' @seealso [getPosteriorFit()], [getRhythmStats()], [getExpectedMeas()]
#'
#' @export
getPosteriorSamples = function(fit, nPosteriorSamples = 200L, overwrite = FALSE) {

  mf = md = NULL
  assertClass(fit, 'limorhyde2')
  assertNumber(nPosteriorSamples, lower = 10)
  nPostSamps = assertCount(nPosteriorSamples, coerce = TRUE)
  assertFlag(overwrite)
  assertTRUE(overwrite || is.null(fit$mashPosteriorSamples))

  mp = foreach(mf = fit$mashFits, md = fit$mashData) %do% {
    mpr = mashr::mash_compute_posterior_matrices(
      mf, md, algorithm.version = 'R', posterior_samples = nPostSamps)
    mpi = mpr$PosteriorSamples}

  co = fit$coefficients
  coArray = array(rep(co, nPostSamps), dim = c(dim(co), nPostSamps),
                  dimnames = dimnames(co))
  coArray[, fit$mashIdx, ] = abind::abind(mp, along = 2L)

  fit$mashPosteriorSamples = coArray
  return(fit)}
