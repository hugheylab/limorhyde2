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
#' * `mashData`: `mash` data object
#' * `mashFit`: `mash` fit object
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

  assertClass(fit, 'limorhyde2')
  covMethod = match.arg(covMethod)
  assertList(getSigResArgs)
  assertCount(npc, positive = TRUE)
  assertList(covEdArgs)
  assertLogical(overwrite, len = 1L)
  assertTRUE(overwrite || is.null(fit$mashFit))

  mashCondCoefs = TRUE
  co = fit$coefficients
  c(shifts, nKnots, nConds) %<-% fit[c('shifts', 'nKnots', 'nConds')]

  idxStart = if (isTRUE(mashCondCoefs)) 2 else nConds + 1
  idxEnd = nConds * (nKnots + 1)
  idxTmp = idxStart:idxEnd # only shrink these

  idx = rep(idxTmp, length(shifts)) +
    rep((0:(length(shifts) - 1)) * ncol(co) / length(shifts),
        each = length(idxTmp))

  md = mashr::mash_set_data(co[, idx], fit$stdErrors[, idx])

  uc = if (covMethod == 'data-driven') NULL else mashr::cov_canonical(md)

  if (covMethod == 'canonical') {
    ued = NULL
  } else {
    m1by1 = mashr::mash_1by1(md)
    strong = do.call(mashr::get_significant_results, c(list(m1by1), getSigResArgs))
    upca = mashr::cov_pca(md, npc = npc, subset = strong)
    ued = do.call(mashr::cov_ed, c(list(md, upca, strong), covEdArgs))}

  mc = mashr::mash(md, c(uc, ued), ...)
  co[, idx] = ashr::get_pm(mc)

  fit$mashData = md
  fit$mashFit = mc
  fit$mashCoefficients = co
  fit$mashIdx = idx
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
getPosteriorSamples = function(fit, nPosteriorSamples = 200, overwrite = FALSE) {

  assertClass(fit, 'limorhyde2')
  assertNumber(nPosteriorSamples, lower = 10)
  nPostSamps = assertCount(nPosteriorSamples, coerce = TRUE)
  assertLogical(overwrite, len = 1L)
  assertTRUE(overwrite || is.null(fit$mashPosteriorSamples))

  mp = mashr::mash_compute_posterior_matrices(
    fit$mashFit, fit$mashData, algorithm.version = 'R',
    posterior_samples = nPostSamps)

  co = fit$coefficients
  coArray = array(rep(co, nPostSamps), dim = c(dim(co), nPostSamps),
                  dimnames = dimnames(co))
  coArray[, fit$mashIdx, ] = mp$PosteriorSamples

  fit$mashPosteriorSamples = coArray
  return(fit)}
