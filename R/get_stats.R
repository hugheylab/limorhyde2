#' @export
getRhythmStats = function(coefMat) {
  c(period, condLevels, nKnots, nConds) %<-%
    attributes(coefMat)[c('period', 'condLevels', 'nKnots', 'nConds')]

  g = function(time) getBasis(time, period, nKnots, TRUE)
  tr = seq(0, period, length.out = nKnots * 20)

  if (nConds == 1L && is.null(condLevels)) {
    condLevels = 'lava'
  } else if (!(nConds > 1 && nConds == length(condLevels))) {
    stop('Inconsistent coefficient matrix and condition levels.')}

  rhyStats = foreach(condIdx = 1:nConds, .combine = rbind) %do% {
    coefNow = getCoefMatOneCond(coefMat, condIdx, nConds, nKnots)

    feo = foreach(co = iterators::iter(coefNow, by = 'row'), .combine = rbind)

    rhyStatsNow = feo %dopar% {
      f = function(time) g(time) %*% t(co)
      d = getOptima(f, tr)
      d[, peak_trough_amp := peak_value - trough_value]
      d[, rms_amp := getRmsAmp(f, co, period)]}

    rhyStatsNow[, mean_value := coefNow[, 1L]]
    rhyStatsNow[, cond := condLevels[condIdx]]
    rhyStatsNow[, feature := rownames(coefMat)]}

  data.table::setcolorder(rhyStats, c('cond', 'feature'))
  rhyStats[, cond := factor(cond, condLevels)]
  if (nConds == 1L) rhyStats[, cond := NULL]
  return(rhyStats[])}


#' @export
getDiffRhythmStats = function(coefMat, rhyStats, condLevels) {
  stopifnot('cond' %in% colnames(rhyStats),
            length(condLevels) == 2L,
            condLevels %in% unique(rhyStats$cond))

  d0 = rhyStats[cond %in% condLevels]
  d0[, cond := factor(cond, condLevels)]
  data.table::setorderv(d0, 'cond')
  period = attr(coefMat, 'period')

  cols = c('mean_value', 'peak_trough_amp', 'rms_amp', 'peak_phase', 'trough_phase')
  diffRhyStats = d0[, lapply(.SD, diff), by = feature, .SDcols = cols]
  diffRhyStats[, peak_phase := centerCircDiff(peak_phase, period)]
  diffRhyStats[, trough_phase := centerCircDiff(trough_phase, period)]
  data.table::setnames(diffRhyStats, cols, paste0('diff_', cols))

  # calculate rms difference in rhythmic fit between conditions
  diffRhyStats[, rms_diff_rhy := getRmsDiffRhy(coefMat, condLevels)]
  return(diffRhyStats[])}


#' @export
getTopTable = function(fit, contrast = c('rhy', 'diff_rhy'), ...) {
  contrast = match.arg(contrast)
  c(nKnots, nConds, nCovars) %<-%
    attributes(fit$coefficients)[c('nKnots', 'nConds', 'nCovars')]
  stopifnot(nConds > 1 || contrast == 'rhy')

  i = if (contrast == 'rhy') 0 else nKnots
  d = limma::topTable(
    fit, coef = (1 + nConds + i):(ncol(fit) - nCovars), number = Inf, ...)
  data.table::setDT(d, keep.rownames = 'feature')
  d = d[, .(feature, pval_rhy = P.Value, qval_rhy = adj.P.Val)]

  if (contrast != 'rhy') {
    data.table::setnames(d, 2:3, c('pval_diff_rhy', 'qval_diff_rhy'))
    d1 = limma::topTable(fit, coef = 2:nConds, number = Inf, ...)
    d1 = data.table::setDT(d1, keep.rownames = 'feature')[
      , .(feature, pval_diff_mean = P.Value, qval_diff_mean = adj.P.Val)]
    d = merge(d, d1, by = 'feature', sort = FALSE)}
  return(d)}
