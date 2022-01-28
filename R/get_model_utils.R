#' @import checkmate
#' @importFrom data.table data.table := set setattr
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom zeallot %<-%
NULL


getShifts = function(nShifts, nKnots, period) {
  if (nKnots == 2L) {
    shifts = 0 # cosinor is invariant to shifts
  } else {
    knotInterval = period / (nKnots + 1)
    shiftInterval = knotInterval / nShifts
    shifts = seq(0, knotInterval - shiftInterval, shiftInterval)}
  return(shifts)}


getMetadata = function(metadata, timeColname, condColname, covarColnames) {
  m = data.table(time = metadata[[timeColname]])

  if (!is.null(condColname)) {
    set(m, j = 'cond', value = factor(metadata[[condColname]]))}

  if (!is.null(covarColnames)) {
    covarColnames = unique(covarColnames)
    for (covarName in covarColnames) {
      set(m, j = paste0('covar_', covarName), value = metadata[[covarName]])}}

  return(m)}


addIntercept = function(b, intercept) {
  if (isTRUE(intercept)) {
    b = cbind(1, b)
    colnames(b)[1L] = 'intercept'}
  return(b)}


getBasis = function(time, period, nKnots, intercept) {
  if (nKnots == 2L) {
    nKnots = 2L
    tt = time / period * 2 * pi
    b = cbind(cos(tt), sin(tt))
  } else {
    knots = seq(0, period, length = nKnots + 2) # including boundary knots
    b = pbs::pbs(time %% period, knots = knots[-c(1, length(knots))],
                 Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
    # scale basis so intercept doesn't change when other coefs shrink
    b = b - 1 / (nKnots + 1)}

  colnames(b) = paste0('basis', 1:nKnots)
  b = addIntercept(b, intercept)
  return(b)}


getDesign = function(metadata, period, nKnots) {
  condIdx = NULL
  b = getBasis(metadata$time, period, nKnots, FALSE)
  m = cbind(metadata, b)

  r = paste0('basis', 1:nKnots, collapse = ' + ')
  if ('cond' %in% colnames(m)) r = sprintf('cond + cond : (%s)', r)

  covarIdx = startsWith(colnames(m), 'covar')
  if (any(covarIdx)) r = paste(c(r, colnames(m)[covarIdx]), collapse = ' + ')
  design = stats::model.matrix(stats::formula(paste('~', r)), data = m)

  covarIdx = startsWith(colnames(design), 'covar')
  if (any(covarIdx)) design = design[, c(which(!covarIdx), which(covarIdx))]

  nConds = length(unique(m$cond)) # works even if cond is not a column
  if (nConds > 1) {
    nCols = ncol(design) - sum(covarIdx)
    idx = foreach(condIdx = 1:nConds, .combine = c) %do% {
      seq(nConds + condIdx, nCols, nConds)}
    idx = c(1:nConds, idx)
    if (any(covarIdx)) idx = c(idx, (nCols + 1):ncol(design))
    design = design[, idx]}
  return(design)}


getNumKnotCondCovar = function(cols) {
  nCovs = sum(startsWith(cols, 'covar'))
  # works for limma and deseq2
  nConds = which(grepl('^cond.+(:|\\.)basis1$', cols))[1L] - 1L
  if (is.na(nConds)) nConds = 1L
  nKnots = as.integer((length(cols) - nConds - nCovs) / nConds)
  return(c(nKnots, nConds, nCovs))}
