# #' @export
# getRhythmStats = function(mat){ # just takes a mat with effects
#
#   c(period, nKnots, cond) %<-% attributes(mat)[-2:-1]
#
#   t = seq(0, period, by = period/(nKnots *10))
#
#   statsAll = foreach(cNum = 1:length(cond), .combine = rbind) %do% {
#
#     cMat = getCond(mat, cNum, length(cond), nKnots)
#     # isolate matrix for specific condition
#     # maybe create function that isolates individual condition matrix by index in original matrix
#
#     statsNow = foreach(mNow = iter(cMat, by = "row"), .combine = rbind) %dopar% {
#
#       funcR = function(x, co = mNow, p = period, nk = nKnots){
#         co2 = matrix(co, ncol = 1)
#         b = getBasis(x, p, nk, intercept = TRUE)
#
#         y = b %*% co2
#
#
#         return(y)
#
#       }
#
#
#       res = getOptimize(funcR, t)
#
#       res[, mean_value := mNow[1]]
#       res[, ampl := peak_value - trough_value]
#       res[, feature := rownames(mNow)]
#
#       # res = data.table(funcR(t))
#
#       return(res) }
#
#     if(length(cond) > 1) {
#
#       statsNow[,cond := cond[cNum]]
#     }
#
#
#     return(statsNow) }
#
#   attr(statsAll, 'period') = period
#   attr(statsAll, 'nKnots') = nKnots
#   attr(statsAll, 'cond') = cond
#
#   return(statsAll)
#
# }

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


# #' @export
# getDiffRhythmStats = function(dat, condIds){
#
#   c(period, nKnots, cond) %<-% tail(attributes(dat),3)
#
#   stopifnot(length(condIds) >=2, condIds %in% dat[, unique(cond)])
#
#   d = dat[cond %in% condIds]
#   d = d[order(factor(cond, levels = condIds))] #preserve order
#
#   cols = setdiff(colnames(d), c('cond', 'feature'))
#   d1 = d[, lapply(.SD, diff), .SDcols = cols, by  = feature]
#   d1[, cond := paste(condIds, collapse = ':')]
#
#   d1[, peak_time := fixDiffPhase(peak_time, period)]
#   d1[, trough_time := fixDiffPhase(trough_time, period)]
#
#   setnames(d1, cols, paste0("diff_", cols))
#
#
#   return(d1)
#
# }

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
