#' @export
getModelFit = function(
  y, metadata, period = 24, nKnots = 4, timeColname = 'time',
  condColname = NULL, covarColnames = NULL,
  method = c('trend', 'voom'), lmFitArgs = list(),
  eBayesArgs = if (method == 'trend') list(trend = TRUE) else list()) {

  stopifnot(ncol(y) == nrow(metadata))
  method = match.arg(method)

  m = getMetadata(metadata, timeColname, condColname, covarColnames)
  design = getDesign(m, period, nKnots)

  v = if (method == 'voom') limma::voom(y, design) else y
  fit = do.call(limma::lmFit, c(list(v, design), lmFitArgs))
  fit = do.call(limma::eBayes, c(list(fit), eBayesArgs))

  # add attributes to coef matrix so the mashed coef matrix has them too
  attr(fit$coefficients, 'period') = period
  attr(fit$coefficients, 'condLevels') = levels(m$cond) # always works

  c(nK, nCon, nCov) %<-% getNumKnotCondCovar(colnames(fit))
  attr(fit$coefficients, 'nKnots') = nK
  attr(fit$coefficients, 'nConds') = nCon
  attr(fit$coefficients, 'nCovars') = nCov
  return(fit)}


#' @export
getMashedCoefs = function(
  fit, mashCondCoefs = TRUE, covMethod = c('data-driven', 'canonical'),
  getSigResArgs = list(), npc = attr(fit$coefficients, 'nKnots'),
  covEdArgs = list(), ...) {

  covMethod = match.arg(covMethod)
  stopifnot(length(mashCondCoefs) == 1L,
            is.logical(mashCondCoefs),
            length(npc) == 1L,
            is.numeric(npc))

  co = fit$coefficients # attributes preserved
  se = sqrt(fit$s2.post) * fit$stdev.unscaled
  c(nKnots, nConds) %<-% attributes(co)[c('nKnots', 'nConds')]

  idxStart = if (isTRUE(mashCondCoefs)) 2 else nConds + 1
  idxEnd = nConds * (nKnots + 1)
  idx = idxStart:idxEnd # only shrink these
  # idx = nConds + 1:(nKnots * nConds) # only shrink these

  md = mashr::mash_set_data(co[, idx], se[, idx])
  uc = mashr::cov_canonical(md)

  if (covMethod == 'data-driven') {
    m1by1 = mashr::mash_1by1(md)
    strong = do.call(mashr::get_significant_results, c(list(m1by1), getSigResArgs))
    upca = mashr::cov_pca(md, npc = npc, subset = strong)
    ued = do.call(mashr::cov_ed, c(list(md, upca, strong), covEdArgs))
  } else {
    ued = NULL}

  mc = mashr::mash(md, c(uc, ued), ...)
  co[, idx] = ashr::get_pm(mc)
  return(co)}
