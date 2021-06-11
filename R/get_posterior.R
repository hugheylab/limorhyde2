#' @export
getPosteriorFit = function(
  fit, covMethod = c('data-driven', 'canonical', 'both'), getSigResArgs = list(),
  npc = fit$nKnots, covEdArgs = list(), overwrite = FALSE, ...) {

  stopifnot(inherits(fit, 'limorhyde2'),
            length(mashCondCoefs) == 1L,
            is.logical(mashCondCoefs),
            length(npc) == 1L,
            is.numeric(npc),
            isTRUE(overwrite) || is.null(fit$mashFit))

  mashCondCoefs = TRUE
  covMethod = match.arg(covMethod)

  co = fit$coefficients
  se = do.call(
    cbind, lapply(fit$lmFits, function(f) sqrt(f$s2.post) * f$stdev.unscaled))
  c(shifts, nKnots, nConds) %<-% fit[c('shifts', 'nKnots', 'nConds')]

  idxStart = if (isTRUE(mashCondCoefs)) 2 else nConds + 1
  idxEnd = nConds * (nKnots + 1)
  idxTmp = idxStart:idxEnd # only shrink these

  idx = rep(idxTmp, length(shifts)) +
    rep((0:(length(shifts) - 1)) * ncol(co) / length(shifts),
        each = length(idxTmp))

  md = mashr::mash_set_data(co[, idx], se[, idx])

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
