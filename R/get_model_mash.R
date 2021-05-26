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
    mNow = data.table::copy(m)
    set(mNow, j = 'time', value = mNow$time + shift)
    design = getDesign(mNow, period, nKnots)

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
  return(fit)}


#' @export
getMashFit = function(
  fit, covMethod = c('data-driven', 'canonical', 'both'), getSigResArgs = list(),
  npc = fit$nKnots, covEdArgs = list(), overwrite = FALSE, ...) {

  mashCondCoefs = TRUE
  covMethod = match.arg(covMethod)

  stopifnot(length(mashCondCoefs) == 1L,
            is.logical(mashCondCoefs),
            length(npc) == 1L,
            is.numeric(npc),
            isTRUE(overwrite) || is.null(fit$mashFit))

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
