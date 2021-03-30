# #' @export
# getModelFit = function(x, metadata, period = 24, timeColname, conditionsColname, nKnots,...){
#
#   sm = getSm(metadata, timeColname, conditionsColname)
#   sm[, cond := as.factor(cond)]
#   bMat = getBasis(sm$time, period, nKnots)
#
#   if(is.null(conditionsColname)){
#     bMat = as.data.table(bMat)
#     formFull = ~ .
#
#     } else{
#       bMat = cbind(sm[, .(cond)], bMat)
#       formFull = ~ cond * .}
#
#   design = stats::model.matrix(formFull, data = bMat)
#
#   fit = limma::lmFit(x, design)
#   fit = limma::eBayes(fit, trend = TRUE,...)
#
#   attr(fit$coefficients, 'period') = period
#   attr(fit$coefficients, 'nKnots') = nKnots
#   attr(fit$coefficients, 'cond') = unique(sm$cond)
#
#
#   return(fit)
# }

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


# #' @export
# getRhythmAsh = function(fit, covMethod = c('canonical', 'data-driven')
#   , getSigResArgs = list(), covEdArgs = list()
#   , npc = attr(fit$coefficients, 'nKnots'), ...){
#
#   c(period, nKnots, cond) %<-% attributes(fit$coefficients)[-2:-1]
#
#   bMat = fit$coefficients
#   idxRemove = getCK(bMat)[1]
#
#   bHat = fit$coefficients[, -(1:idxRemove)]
#   sHat = sqrt(fit$s2.post) * fit$stdev.unscaled[, -(1:idxRemove)]
#
#   data = mashr::mash_set_data(bHat, sHat)
#   Uc = mashr::cov_canonical(data)
#
#   covType = match.arg(covMethod)
#
#   if ('data-driven' %in% covMethod
#     & !is.null(npc)) {
#
#     m1by1 = mashr::mash_1by1(data)
#
#     strong = do.call(mashr::get_significant_results
#       , c(list(m1by1), getSigResArgs))
#
#     Upca = mashr::cov_pca(data = data, subset = strong, npc = npc)
#
#     Ued = do.call(mashr::cov_ed, c(list(data = data), list(Ulist_init = Upca)
#       , list(subset = strong), covEdArgs))
#
#   } else if ('data-driven' %in% covMethod
#       & is.null(npc)) {
#
#     stop("Data-driven method specified without argument 'npc'.")
#
#   } else { Ued = NULL }
#
#   resMash = mashr::mash(data,c(Uc, Ued))
#   pm = resMash$result$PosteriorMean
#
#   pm = cbind(bMat[, 1:idxRemove, drop = FALSE], pm)
#
#   attr(pm, 'period') = period
#   attr(pm, 'nKnots') = nKnots
#   attr(pm, 'cond') = cond
#
#   return(pm)}

#' @export
getMashedCoefs = function(
  fit, mashCondCoefs = FALSE, covMethod = c('data-driven', 'canonical'),
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
