#' Calculate rhythmic statistics from a fitted model
#'
#' `getRhythmStats`returns the rhythmic statistics calculated under a `limorhyde2` model for selected features.
#'
#' @param fit A limorhyde2 fit object, as provided by `getModelFit` or `getPosteriorFit`.
#' @param fitType String with potential values of 'posterior_mean', 'posterior_samples', or 'raw' indicating the type of `limorhyde2` object on which to calculate rhythmic statistics.
#' @param features Vector containing the names, row numbers, or logical conditions for features on which to calculate rhythmic statistics.
#'
#' @return A data.table with rows for each feature and columns for the following rhythmic statistics
#' * mean value: mean y value for given feature
#' * peak-to-trough amplitude: peak amplitude - trough amplitude
#' * root mean-squared amplitude: square root of mean squared difference between fitted time and intercept integrated across period length for given feature
#' * peak phase: phase at which peak y-value occurs
#' * trough phase: phase at which minimum y-value occurs
#' as well as condition and posterior sample, if applicable.
#'
#' @export
getRhythmStats = function(
  fit, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL) {

  stopifnot(inherits(fit, 'limorhyde2'))

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


#' Calculate differential rhythmic statistics between conditions
#'
#' `getDiffRhythmStats` returns differences in rhythmic statistics calculated under pairs of
#'  conditions.
#'
#' @param fit A limorhyde2 fit object, as provided by `getModelFit` or `getPosteriorFit`.
#' @param rhyStats A data.table of rhythmic statistics, as returned by `getRhythmStats`. The data.table returned by `getDiffRhythmStats` inherits its `fitType` from `rhyStats`.
#' @param condLevels A character vector containing the two conditions to compare.
#'
#' @return A data.table of the following differential rhythmic statistics
#' * differential mean value: difference in mean value of y-value between features under paired conditions
#' * differential peak-trough amplitude: difference in peak-trough amplitude between paired conditions
#' * differential root mean-squared amplitude: square root of mean squared difference between y-values in paired conditions integrated across period length
#' * differential peak phase: circular difference in peak phase between paired conditions
#' * differential trough phase: circular difference in trough phase between paired conditions.
#'
#' @export
getDiffRhythmStats = function(fit, rhyStats, condLevels) {
  stopifnot(inherits(fit, 'limorhyde2'),
            isTRUE(attr(rhyStats, 'statType') == 'rhy'),
            'cond' %in% colnames(rhyStats),
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
  attr(diffRhyStats, 'condLevels') = condLevels
  return(diffRhyStats[])}


#' Calculate credible intervals for rhythmic statistics
#'
#' `getStatsIntervals` constructs credible intervals for each of the rhythmic statistics calculated in `getRhythmStats` or differential rhythmic statistics calculated in `getDiffRhythmStats`.
#'
#' @param posteriorStats A data.table of posterior samples of (differential) rhythmic statistics from `getRhythmStats` or `getDiffRhythmStats`.
#' @param mass The probability mass for which to calculate the interval.
#' @param method String for type of interval: 'eti' for equal-tailed or 'hdi' for highest (posterior) density. Equal-tailed intervals use [stats::quantiles], while HPD intervals use [HDInterval::hdi].
#'
#' @return A data.table with rows for each feature and (differential) rhythmic statistic, as well as condition if applicable. Columns correspond to the upper and lower bounds of the credible interval.
#'
#' @seealso [getRhythmStats], [getDiffRhythmStats], [stats::quantiles] ,[HDInterval::hdi]
#'
#' @export
getStatsIntervals = function(
  posteriorStats, mass = 0.9, method = c('eti', 'hdi')) {
  # TODO: extend for phase-based stats, possibly in 2D

  stopifnot(isTRUE(attr(posteriorStats, 'fitType') == 'posterior_samples'),
            isTRUE(attr(posteriorStats, 'statType') %in% c('rhy', 'diff_rhy')))
  method = match.arg(method)

  statType = attr(posteriorStats, 'statType')
  idCols = c('feature', 'posterior_sample')
  if (statType == 'diff_rhy') idCols = c('cond', idCols)

  varName = 'statistic'
  byCols = c(setdiff(idCols, 'posterior_sample'), varName)

  measCols = if (statType == 'rhy') {
    c('peak_value', 'trough_value', 'peak_trough_amp', 'rms_amp')
  } else {
    c('diff_mean_value', 'diff_peak_trough_amp', 'diff_rms_amp', 'rms_diff_rhy')}

  d1 = data.table::melt(
    posteriorStats, id.vars = idCols, measure.vars = measCols,
    variable.name = varName, variable.factor = FALSE)
  getInterval = if (method == 'eti') getEti else getHdi
  d2 = d1[, getInterval(value, mass), by = byCols]

  attr(d2, 'statType') = statType
  attr(d2, 'mass') = mass
  attr(d2, 'method') = method
  return(d2)}
