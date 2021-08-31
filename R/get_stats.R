#' Compute rhythmic statistics from fitted models
#'
#' This function uses [stats::optim()] to compute various properties of
#' fitted curves with respect to time, potentially in each condition and for
#' each posterior sample, and adjusting for any covariates.
#'
#' @param fit A `limorhyde2` object.
#' @param fitType String indicating which fitted models to use to compute the
#'   rhythmic statistics. A typical analysis using `limorhyde2` will be based on
#'   'posterior_mean', the default.
#' @param features Vector of names, row numbers, or logical values for
#'   subsetting the features. `NULL` indicates all features.
#' @param dopar Logical indicating whether to run calculations in parallel if
#'   a parallel backend is already set up, e.g., using
#'   [doParallel::registerDoParallel()]. Recommended to minimize runtime.
#'
#' @return A `data.table` containing the following rhythmic statistics:
#'
#' * `peak_phase`: time between 0 and `fit$period` at which the peak or maximum
#'   value occurs
#' * `peak_value`
#' * `trough_phase`: time between 0 and `fit$period` at which the trough or
#'   minimum value occurs
#' * `trough_value`
#' * `peak_trough_amp`: `peak_value - trough_value`
#' * `rms_amp`: root mean square difference between fitted curve and mean value
#'   between time 0 and `fit$period`
#' * `mean_value`: between time 0 and `fit$period`
#'
#' The rows of the `data.table` depend on the `fit` object and `fitType`:
#'
#' * `fit` contains data from one condition and `fitType` is posterior_mean' or
#'   'raw': one row per feature.
#' * `fit` contains data from one condition and `fitType` is
#'   'posterior_samples': one row per feature per posterior sample.
#'  * `fit` contains data from multiple conditions and `fitType` is
#'   'posterior_mean' or 'raw': one row per feature per condition.
#'  * `fit` contains data from multiple conditions and `fitType` is
#'   'posterior_samples': one row per feature per condition per posterior
#'   sample.
#'
#' @seealso [getModelFit()], [getPosteriorFit()],  [getPosteriorSamples()],
#'   [getDiffRhythmStats()], [getStatsIntervals()]
#'
#' @export
getRhythmStats = function(
  fit, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL, dopar = TRUE) {

  shifts = period = nKnots = nConds = postSampIdx = condIdx = cond = NULL
  peak_value = trough_value = peak_trough_amp = co = posterior_sample = NULL

  assertClass(fit, 'limorhyde2')
  fitType = match.arg(fitType)
  checkFitType(fit, fitType)
  assertLogical(dopar, any.missing = FALSE, len = 1L)

  c(shifts, period, conds, nKnots, nConds) %<-%
    fit[c('shifts', 'period', 'conds', 'nKnots', 'nConds')]

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, TRUE)}))}

  tr = seq(0, period, length.out = nKnots * 20)
  if (nConds == 1L) conds = 'lava'

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  doPost = if (nPostSamps == 1L | !dopar) `%do%` else `%dopar%`
  doFeat = if (nPostSamps == 1L & dopar) `%dopar%` else `%do%`

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
        set(d, j = 'rms_amp', value = getRmsAmp(f, co, period))})

      idx = seq(1, ncol(coefNow), ncol(coefNow) / length(shifts))
      set(r2, j = 'mean_value', value = rowMeans(coefNow[, idx, drop = FALSE]))
      set(r2, j = 'cond', value = conds[condIdx])
      set(r2, j = 'feature', value = rownames(coefMat))}

    set(r1, j = 'posterior_sample', value = postSampIdx)})

  data.table::setcolorder(rhyStats, c('cond', 'feature', 'posterior_sample'))
  rhyStats[, cond := factor(cond, conds)]
  if (nConds == 1L) rhyStats[, cond := NULL]
  if (nPostSamps == 1L) rhyStats[, posterior_sample := NULL]

  setattr(rhyStats, 'statType', 'rhy')
  setattr(rhyStats, 'fitType', fitType)
  return(rhyStats)}


#' Compute differentially rhythmic statistics from fitted models
#'
#' This function computes differences in rhythmicity between fitted curves for a
#' given pair of conditions.
#'
#' @param fit A `limorhyde2` object containing data from multiple conditions.
#' @param rhyStats A `data.table` of rhythmic statistics, as returned by
#'   [getRhythmStats()], for fitted models in `fit`.
#' @param conds A character vector indicating the conditions to compare
#'   pairwise, by default all conditions in `fit`.
#' @param dopar Logical indicating whether to run calculations in parallel if
#'   a parallel backend is already set up, e.g., using
#'   [doParallel::registerDoParallel()]. Recommended to minimize runtime.
#'
#' @return A `data.table` containing the following differentially rhythmic
#'   statistics:
#'
#' * `diff_mean_value`
#' * `diff_peak_trough_amp`
#' * `diff_rms_amp`
#' * `diff_peak_phase`: circular difference between `-fit$period/2` and
#'   `fit$period/2`
#' * `diff_trough_phase`: circular difference between `-fit$period/2` and
#'   `fit$period/2`
#' * `rms_diff_rhy`: root mean square difference in mean-centered fitted curves
#'
#' The stats will be based on the value for `cond2` minus the value for `cond1`.
#' The rows of the `data.table` depend on the 'fitType' attribute of `rhyStats`:
#'
#' * 'fitType' is 'posterior_mean' or 'raw': one row per feature per pair of
#'   conditions.
#' * 'fitType' is 'posterior_samples': one row per feature per posterior sample
#'   per pair of conditions.
#'
#' @seealso [getRhythmStats()], [getStatsIntervals()]
#'
#' @export
getDiffRhythmStats = function(fit, rhyStats, conds = fit$conds, dopar = TRUE) {

  cond = .SD = diff_peak_phase = diff_trough_phase = condInt = . = cond1  =
    condInt1 = mean_value1 = peak_trough_amp1 = rms_amp1 = peak_phase1 =
    trough_phase1 = cond2 = condInt2 = mean_value2 = peak_trough_amp2 =
    rms_amp2 = peak_phase2 = trough_phase2 = NULL

  conds = unique(conds)
  assertClass(fit, 'limorhyde2')
  assertTRUE(fit$nConds >= 2)
  assertDataTable(rhyStats)
  assertTRUE(attr(rhyStats, 'statType') == 'rhy')
  assertTRUE('cond' %in% colnames(rhyStats))
  assertSubset(conds, fit$conds)
  assertSubset(conds, levels(rhyStats$cond))
  assertLogical(dopar, any.missing = FALSE, len = 1L)

  d0 = rhyStats[cond %in% conds]
  set(d0, j = 'cond', value = factor(d0$cond, conds))
  d0[, condInt := as.integer(cond)]

  fitType = attr(rhyStats, 'fitType')
  byCols = c('feature', if (fitType == 'posterior_samples') 'posterior_sample')

  diffRhyStats = merge(
    d0, d0, allow.cartesian = TRUE, by = byCols, suffixes = 1:2)
  diffRhyStats = diffRhyStats[condInt1 < condInt2]
  diffRhyStats = diffRhyStats[,
    .(cond1, cond2,
      diff_mean_value = mean_value2 - mean_value1,
      diff_peak_trough_amp = peak_trough_amp2 - peak_trough_amp1,
      diff_rms_amp = rms_amp2 - rms_amp1,
      diff_peak_phase = peak_phase2 - peak_phase1,
      diff_trough_phase = trough_phase2 - trough_phase1),
    by = byCols]

  diffRhyStats[, diff_peak_phase := centerCircDiff(diff_peak_phase, fit$period)]
  diffRhyStats[, diff_trough_phase := centerCircDiff(diff_trough_phase, fit$period)]

  # calculate rms difference in rhythmic fit between conditions
  pairs = unique(diffRhyStats[, .(cond1, cond2)])
  featureIdx = rownames(fit$coefficients) %in% unique(rhyStats$feature)

  feo = foreach(cond1 = pairs$cond1, cond2 = pairs$cond2, .combine = rbind)
  rmsDiffRhy = feo %do% {
    rmsDiffRhyTmp = getRmsDiffRhy(fit, c(cond1, cond2), fitType, featureIdx, dopar)
    rmsDiffRhyTmp = data.table(rmsDiffRhyTmp, cond1, cond2)}

  diffRhyStats = merge(
    diffRhyStats, rmsDiffRhy, sort = FALSE, by = c(byCols, 'cond1', 'cond2'))

  setattr(diffRhyStats, 'statType', 'diff_rhy')
  setattr(diffRhyStats, 'fitType', fitType)
  return(diffRhyStats)}


#' Compute credible intervals for rhythmic or differentially rhythmic statistics
#'
#' This function uses posterior samples to quantify uncertainty in the
#' properties of fitted curves.
#'
#' @param posteriorStats A `data.table` of statistics for posterior samples, as
#'   returned by [getRhythmStats()] or [getDiffRhythmStats()].
#' @param mass Number between 0 and 1 indicating the probability mass for which
#'   to calculate the intervals.
#' @param method String indicating the type of interval: 'eti' for equal-tailed
#'   using [stats::quantile()], or 'hdi' for highest density using
#'   [HDInterval::hdi()].
#'
#' @return A `data.table` containing lower and upper bounds of various
#'   statistics for each feature or each feature-condition pair.
#'
#' @seealso [getRhythmStats()], [getDiffRhythmStats()],
#'   [getExpectedMeasIntervals()]
#'
#' @export
getStatsIntervals = function(
  posteriorStats, mass = 0.9, method = c('eti', 'hdi')) {
  # TODO: extend for phase-based stats, possibly in 2D

  value = NULL
  assertDataTable(posteriorStats)
  assertTRUE(attr(posteriorStats, 'fitType') == 'posterior_samples')
  assertChoice(attr(posteriorStats, 'statType'), c('rhy', 'diff_rhy'))
  assertNumber(mass, lower = 0.5, upper = 1 - .Machine$double.eps)
  method = match.arg(method)

  statType = attr(posteriorStats, 'statType')
  condCols = if (statType == 'rhy') 'cond' else c('cond1', 'cond2')
  idCols = intersect(c(condCols, 'feature', 'posterior_sample'),
                     colnames(posteriorStats))

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

  setattr(d2, 'statType', statType)
  setattr(d2, 'mass', mass)
  setattr(d2, 'method', method)
  return(d2)}
