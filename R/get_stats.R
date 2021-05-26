#' @export
getRhythmStats = function(
  fit, coefType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  c(shifts, period, condLevels, nKnots, nConds) %<-%
    fit[c('shifts', 'period', 'condLevels', 'nKnots', 'nConds')]

  coefType = match.arg(coefType)
  if (coefType == 'posterior_mean' && is.null(fit$mashCoefficients)) {
    stop('No mash coefficients from which to calculate statistics.')
  } else if (coefType == 'posterior_samples' && is.null(fit$mashPosteriorSamples)) {
    stop('No mash posterior samples from which to calculate statistics.')}

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, TRUE)}))}

  tr = seq(0, period, length.out = nKnots * 20)

  if (nConds == 1L && is.null(condLevels)) {
    condLevels = 'lava'
  } else if (!(nConds > 1 && nConds == length(condLevels))) {
    stop('Inconsistent coefficient matrix and condition levels.')}

  coefArray = getCoefArray(fit, coefType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  doPost = if (nPostSamps == 1L) `%do%` else `%dopar%`
  doFeat = if (nPostSamps == 1L) `%dopar%` else `%do%`

  rhyStats = doPost(foreach(postSampIdx = 1:nPostSamps, .combine = rbind), {
    coefMat = abind::adrop(coefArray[, , postSampIdx, drop = FALSE], drop = 3)

    r1 = foreach(condIdx = 1:nConds, .combine = rbind) %do% {
      coefNow = getCoefMatOneCond(
        coefMat, condIdx, nConds, nKnots, length(shifts))

      feo = foreach(co = iterators::iter(coefNow, by = 'row'), .combine = rbind)
      r2 = doFeat(feo, {
        f = function(time) (g(time) %*% t(co)) / length(shifts)
        d = getOptima(f, tr)
        d[, peak_trough_amp := peak_value - trough_value]
        d[, rms_amp := getRmsAmp(f, co, period)]})

      idx = seq(1, ncol(coefNow), ncol(coefNow) / length(shifts))
      set(r2, j = 'mean_value', value = rowMeans(coefNow[, idx, drop = FALSE]))
      set(r2, j = 'cond', value = condLevels[condIdx])
      set(r2, j = 'feature', value = rownames(coefMat))}

    set(r1, j = 'posterior_sample', value = postSampIdx)})

  data.table::setcolorder(rhyStats, c('cond', 'feature', 'posterior_sample'))
  rhyStats[, cond := factor(cond, condLevels)]
  if (nConds == 1L) rhyStats[, cond := NULL]
  if (nPostSamps == 1L) rhyStats[, posterior_sample := NULL]

  attr(rhyStats, 'statType') = 'rhy'
  attr(rhyStats, 'coefType') = coefType
  return(rhyStats[])}


#' @export
getDiffRhythmStats = function(fit, rhyStats, condLevels) {
  stopifnot('cond' %in% colnames(rhyStats),
            length(condLevels) == 2L,
            all(condLevels %in% fit$condLevels),
            all(condLevels %in% unique(rhyStats$cond)))

  d0 = rhyStats[cond %in% condLevels]
  d0[, cond := factor(cond, condLevels)]
  data.table::setorderv(d0, 'cond')

  coefType = attr(rhyStats, 'coefType')
  byCols = c('feature', if (coefType == 'posterior_samples') 'posterior_sample')
  cols = c('mean_value', 'peak_trough_amp', 'rms_amp', 'peak_phase', 'trough_phase')

  diffRhyStats = d0[, lapply(.SD, diff), by = byCols, .SDcols = cols]
  diffRhyStats[, peak_phase := centerCircDiff(peak_phase, fit$period)]
  diffRhyStats[, trough_phase := centerCircDiff(trough_phase, fit$period)]
  data.table::setnames(diffRhyStats, cols, paste0('diff_', cols))

  # calculate rms difference in rhythmic fit between conditions
  featureIdx = rownames(fit$coefficients) %in% unique(rhyStats$feature)
  rmsDiffRhy = getRmsDiffRhy(fit, condLevels, coefType, featureIdx)
  diffRhyStats = merge(diffRhyStats, rmsDiffRhy, sort = FALSE)

  attr(diffRhyStats, 'statType') = 'diff_rhy'
  attr(diffRhyStats, 'coefType') = coefType
  return(diffRhyStats[])}
