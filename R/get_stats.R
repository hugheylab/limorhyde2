#' @export
getRhythmStats = function(
  fit, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  fitType = match.arg(fitType)
  if (fitType == 'posterior_mean' && is.null(fit$mashCoefficients)) {
    stop('No posterior mean to calculate statistics, please run getPosteriorFit.')
  } else if (fitType == 'posterior_samples' && is.null(fit$mashPosteriorSamples)) {
    stop('No posterior samples to calculate statistics, please run getPosteriorSamples.')}

  c(shifts, period, condLevels, nKnots, nConds) %<-%
    fit[c('shifts', 'period', 'condLevels', 'nKnots', 'nConds')]

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, TRUE)}))}

  tr = seq(0, period, length.out = nKnots * 20)
  if (nConds == 1L) condLevels = 'lava'

  coefArray = getCoefArray(fit, fitType)
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
  attr(rhyStats, 'fitType') = fitType
  return(rhyStats[])}


#' @export
getDiffRhythmStats = function(fit, rhyStats, condLevels) {
  stopifnot('cond' %in% colnames(rhyStats),
            attr(rhyStats, 'statType') == 'rhy',
            length(condLevels) == 2L,
            all(condLevels %in% fit$condLevels),
            all(condLevels %in% unique(rhyStats$cond)))

  d0 = rhyStats[cond %in% condLevels]
  d0[, cond := factor(cond, condLevels)]
  data.table::setorderv(d0, 'cond')

  fitType = attr(rhyStats, 'fitType')
  byCols = c('feature', if (fitType == 'posterior_samples') 'posterior_sample')
  cols = c('mean_value', 'peak_trough_amp', 'rms_amp', 'peak_phase', 'trough_phase')

  diffRhyStats = d0[, lapply(.SD, diff), by = byCols, .SDcols = cols]
  diffRhyStats[, peak_phase := centerCircDiff(peak_phase, fit$period)]
  diffRhyStats[, trough_phase := centerCircDiff(trough_phase, fit$period)]
  data.table::setnames(diffRhyStats, cols, paste0('diff_', cols))

  # calculate rms difference in rhythmic fit between conditions
  featureIdx = rownames(fit$coefficients) %in% unique(rhyStats$feature)
  rmsDiffRhy = getRmsDiffRhy(fit, condLevels, fitType, featureIdx)
  diffRhyStats = merge(diffRhyStats, rmsDiffRhy, sort = FALSE)

  attr(diffRhyStats, 'statType') = 'diff_rhy'
  attr(diffRhyStats, 'fitType') = fitType
  return(diffRhyStats[])}
