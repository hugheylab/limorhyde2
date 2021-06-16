#' Perform multivariate adaptive shrinkage (mash) on a fitted model
#'
#' Given a linear model fit, from `getModelFit`, `getPosteriorFit` fits a mash model to coefficients to improve coefficient estimates. The mash model returns moderated coefficients as Posterior Means.
#'
#' @param fit A `limorhyde2` fit object, as provided by `getModelFit` or `getPosteriorFit`.
#' @param covMethod String indicating types of covariance matrices to compute when fitting the `mash` model. Can take one of 'data-driven', 'canonical', or 'both'
#' @param getSigResArgs List of argument to be passed to [mashr::get_significant_results]. Used to find significant effects from at least one condition.
#' @param npc Number of principle components to use in principal component analysis of fitted model
#' @param covEdArgs List specifying argument-value pairs to be passed to [mashr::cov_ed].
#' @param overwrite Logical for whether to recompute mashr fit if it already exists.
#'
#' @return A `limorhyde2` class object containing everything found in `fit` with added elements:
#'
#' * `mashData` a data object of class `mash`
#' * `mashFit` list of mash fit results
#' * `mashCoefficients` matrix of resulting posterior model coefficients
#' * `mashIdx` vector of fit coefficient indices included in `mash`fitting
#'
#' @seealso [getModelFit], [mashr::get_significant_results], [mashr::cov_ed]
#'
#' @export
getPosteriorFit = function(
  fit, covMethod = c('data-driven', 'canonical', 'both'), getSigResArgs = list(),
  npc = fit$nKnots, covEdArgs = list(), overwrite = FALSE, ...) {

  stopifnot(inherits(fit, 'limorhyde2'),
            length(npc) == 1L,
            is.numeric(npc),
            isTRUE(overwrite) || is.null(fit$mashFit))

  mashCondCoefs = TRUE
  covMethod = match.arg(covMethod)

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


#' Draw samples from posterior estimate distributions from a fitted model
#'
#' After applying a mash model with `getPosteriorFit`, which computes posterior distribution of each effect, use `getPosteriorSamples` to sample from the posterior for a set of estimates for each coefficient.
#'
#' @param fit A limorhyde2 object, as provided by `getModelFit` or `getPosteriorFit`.
#' @param nPosteriorSamples number of samples to draw from the posterior distribution of each effect
#' @param overwrite Logical indicating whether to posterior samples should be recomputed.
#'
#' @return a limorhyde2 object containing everything found in `fit` with added element:
#' * `mashPosteriorSamples`: An m x n x p array with m features, n model coefficients, and p posterior samples
#'
#' @export
getPosteriorSamples = function(fit, nPosteriorSamples = 200, overwrite = FALSE) {

  stopifnot(!is.null(fit$mashFit),
            isTRUE(overwrite) || is.null(fit$mashPosteriorSamples),
            length(nPosteriorSamples) == 1L,
            is.numeric(nPosteriorSamples),
            nPosteriorSamples >= 10)

  nPostSamps = round(nPosteriorSamples)

  mp = mashr::mash_compute_posterior_matrices(
    fit$mashFit, fit$mashData, algorithm.version = 'R',
    posterior_samples = nPostSamps)

  co = fit$coefficients
  coArray = array(rep(co, nPostSamps), dim = c(dim(co), nPostSamps),
                  dimnames = dimnames(co))
  coArray[, fit$mashIdx, ] = mp$PosteriorSamples

  fit$mashPosteriorSamples = coArray
  return(fit)}
