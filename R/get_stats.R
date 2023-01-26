#' Compute rhythm statistics from fitted models
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
#' @param rms Logical indicating whether to calculate `rms_amp`.
#'
#' @return A `data.table` containing the following rhythm statistics:
#'
#' * `peak_phase`: time between 0 and `fit$period` at which the peak or maximum
#'   value occurs
#' * `peak_value`
#' * `trough_phase`: time between 0 and `fit$period` at which the trough or
#'   minimum value occurs
#' * `trough_value`
#' * `peak_trough_amp`: `peak_value - trough_value`
#' * `rms_amp`: root mean square difference between fitted curve and mean value
#'   between time 0 and `fit$period` (only calculated if `rms` is `TRUE`)
#' * `mesor`: mean value between time 0 and `fit$period`
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
#' @eval examples1()
#'
#' @seealso [getModelFit()], [getPosteriorFit()],  [getPosteriorSamples()],
#'   [getDiffRhythmStats()], [getStatsIntervals()]
#'
#' @export
getRhythmStats = function(
  fit, fitType = c('posterior_mean', 'posterior_samples', 'raw'),
  features = NULL, dopar = TRUE, rms = FALSE) {

  shifts = period = nKnots = degree = nConds = postSampIdx = condIdx = cond =
    peak_value = trough_value = peak_trough_amp = co = posterior_sample = NULL

  assertClass(fit, 'limorhyde2')
  fitType = match.arg(fitType)
  checkFitType(fit, fitType)
  assertFlag(dopar)
  assertFlag(rms)

  c(shifts, period, conds, nKnots, degree, nConds) %<-%
    fit[c('shifts', 'period', 'conds', 'nKnots', 'degree', 'nConds')]

  g = function(time) {
    do.call(cbind, lapply(shifts, function(shift) {
      getBasis(time + shift, period, nKnots, degree, TRUE)}))}

  tr = seq(0, period, length.out = nKnots * 20)
  if (nConds == 1L) conds = 'lava'

  coefArray = getCoefArray(fit, fitType)
  nPostSamps = dim(coefArray)[3L]
  if (!is.null(features)) coefArray = coefArray[features, , , drop = FALSE]

  reg = foreach::getDoParRegistered()
  if (dopar && isAlreadyInParallel()) {
    warning('getRhythmStats was called with dopar = TRUE within a separate parallel call. Setting dopar to FALSE.')
    dopar = FALSE}
  doPost = if (nPostSamps == 1L || !dopar || !reg) `%do%` else `%dopar%`
  doFeat = if (nPostSamps == 1L && dopar && reg) `%dopar%` else `%do%`

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
        if (rms) set(d, j = 'rms_amp', value = getRmsAmp(f, co, period))
        d})

      idx = seq(1, ncol(coefNow), ncol(coefNow) / length(shifts))
      set(r2, j = 'mesor', value = rowMeans(coefNow[, idx, drop = FALSE]))
      set(r2, j = 'cond', value = conds[condIdx])
      set(r2, j = 'feature', value = rownames(coefMat))}

    set(r1, j = 'posterior_sample', value = postSampIdx)})

  data.table::setcolorder(rhyStats, c('cond', 'feature', 'posterior_sample'))
  rhyStats[, cond := factor(cond, conds)]
  if (nConds == 1L) rhyStats[, cond := NULL]
  if (nPostSamps == 1L) rhyStats[, posterior_sample := NULL]

  setattr(rhyStats, 'statType', 'rhy')
  setattr(rhyStats, 'fitType', fitType)
  setattr(rhyStats, 'period', period)
  return(rhyStats)}


#' Compute differential rhythm statistics from fitted models
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
#' @return A `data.table` containing the following differential rhythm
#'   statistics:
#'
#' * `mean_mesor`
#' * `mean_peak_trough_amp`
#' * `mean_rms_amp` (only calculated if `rms` to [getRhythmStats()] was `TRUE`)
#' * `diff_mesor`
#' * `diff_peak_trough_amp`
#' * `diff_rms_amp` (only calculated if `rms` to [getRhythmStats()] was `TRUE`)
#' * `diff_peak_phase`: circular difference between `-fit$period/2` and
#'   `fit$period/2`
#' * `diff_trough_phase`: circular difference between `-fit$period/2` and
#'   `fit$period/2`
#' * `diff_rhy_dist`: Euclidean distance between polar coordinates
#'   (`peak_trough_amp`, `peak_phase`)
#' * `rms_diff_rhy`: root mean square difference in mean-centered fitted curves
#'    (only calculated if `rms` to [getRhythmStats()] was `TRUE`)
#'
#' The stats will be based on the value for `cond2` minus the value for `cond1`.
#' The rows of the `data.table` depend on the 'fitType' attribute of `rhyStats`:
#'
#' * 'fitType' is 'posterior_mean' or 'raw': one row per feature per pair of
#'   conditions.
#' * 'fitType' is 'posterior_samples': one row per feature per posterior sample
#'   per pair of conditions.
#'
#' @eval examples1()
#'
#' @seealso [getRhythmStats()], [getStatsIntervals()]
#'
#' @export
getDiffRhythmStats = function(
  fit, rhyStats, conds = fit$conds, dopar = TRUE) {

  cond = .SD = diff_peak_phase = diff_trough_phase = condInt = . = cond1  =
    condInt1 = mesor1 = peak_trough_amp1 = rms_amp1 = peak_phase1 =
    trough_phase1 = cond2 = condInt2 = mesor2 = peak_trough_amp2 =
    rms_amp2 = peak_phase2 = trough_phase2 = mean_rms_amp = diff_rms_amp = NULL

  conds = unique(conds)
  assertClass(fit, 'limorhyde2')
  assertTRUE(fit$nConds >= 2)
  assertDataTable(rhyStats)
  assertTRUE(attr(rhyStats, 'statType') == 'rhy')
  assertTRUE('cond' %in% colnames(rhyStats))
  assertSubset(conds, fit$conds)
  assertSubset(conds, levels(rhyStats$cond))
  assertFlag(dopar)

  d0 = rhyStats[cond %in% conds]
  set(d0, j = 'cond', value = factor(d0$cond, conds))
  d0[, condInt := as.integer(cond)]

  period = fit$period
  fitType = attr(rhyStats, 'fitType')
  byCols = c('feature', if (fitType == 'posterior_samples') 'posterior_sample')

  diffRhyStats0 = merge(
    d0, d0, allow.cartesian = TRUE, by = byCols, suffixes = 1:2)
  diffRhyStats0 = diffRhyStats0[condInt1 < condInt2]

  diffRhyStats = diffRhyStats0[
    , .(cond1, cond2,
        mean_mesor = 0.5 * (mesor1 + mesor2),
        mean_peak_trough_amp = 0.5 * (peak_trough_amp1 + peak_trough_amp2),
        diff_mesor = mesor2 - mesor1,
        diff_peak_trough_amp = peak_trough_amp2 - peak_trough_amp1,
        diff_peak_phase = peak_phase2 - peak_phase1,
        diff_trough_phase = trough_phase2 - trough_phase1,
        diff_rhy_dist = getDist(peak_trough_amp1, peak_phase1,
                                peak_trough_amp2, peak_phase2, period)),
    by = byCols]

  if ('rms_amp' %in% colnames(rhyStats)) {
    diffRhyStatsRms = diffRhyStats0[
      , .(cond1, cond2,
          mean_rms_amp = 0.5 * (rms_amp1 + rms_amp2),
          diff_rms_amp = rms_amp2 - rms_amp1),
      by = byCols]
    diffRhyStats = cbind(
      diffRhyStats, diffRhyStatsRms[, .(mean_rms_amp, diff_rms_amp)])}

  diffRhyStats[, diff_peak_phase := centerCircDiff(diff_peak_phase, period)]
  diffRhyStats[, diff_trough_phase := centerCircDiff(diff_trough_phase, period)]

  # calculate rms difference in rhythmic fit between conditions
  if ('rms_amp' %in% colnames(rhyStats)) {
    pairs = unique(diffRhyStats[, .(cond1, cond2)])
    featureIdx = rownames(fit$coefficients) %in% unique(rhyStats$feature)

    feo = foreach(cond1 = pairs$cond1, cond2 = pairs$cond2, .combine = rbind)
    rmsDiffRhy = feo %do% {
      rmsDiffRhyTmp = getRmsDiffRhy(
        fit, c(as.character(cond1), as.character(cond2)), fitType, featureIdx, dopar)
      rmsDiffRhyTmp = data.table(rmsDiffRhyTmp, cond1, cond2)}

    diffRhyStats = merge(
      diffRhyStats, rmsDiffRhy, sort = FALSE, by = c(byCols, 'cond1', 'cond2'))}

  setattr(diffRhyStats, 'statType', 'diff_rhy')
  setattr(diffRhyStats, 'fitType', fitType)
  setattr(diffRhyStats, 'period', period)
  return(diffRhyStats)}


#' Compute credible intervals for rhythm or differential rhythm statistics
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
#'   statistics for each feature or each feature-condition pair. For
#'   `peak_trough_amp` and `rms_amp`, a negative lower bound indicates a rhythm
#'   of the opposite phase.
#'
#' @eval examples2()
#'
#' @seealso [getRhythmStats()], [getDiffRhythmStats()],
#'   [getExpectedMeasIntervals()]
#'
#' @export
getStatsIntervals = function(
  posteriorStats, mass = 0.9, method = c('eti', 'hdi')) {

  value = . = amp = mean_phase = peak_phase = NULL
  assertDataTable(posteriorStats)
  assertTRUE(attr(posteriorStats, 'fitType') == 'posterior_samples')
  assertChoice(attr(posteriorStats, 'statType'), c('rhy', 'diff_rhy'))
  assertNumber(attr(posteriorStats, 'period'), lower = .Machine$double.eps)
  assertNumber(mass, lower = 0.5, upper = 1 - .Machine$double.eps)
  method = match.arg(method)

  statType = attr(posteriorStats, 'statType')
  period = attr(posteriorStats, 'period')

  condCols = if (statType == 'rhy') 'cond' else c('cond1', 'cond2')
  idCols = intersect(c(condCols, 'feature', 'posterior_sample'),
                     colnames(posteriorStats))

  varName = 'statistic'
  byCols = c(setdiff(idCols, 'posterior_sample'), varName)

  measCols = if (statType == 'rhy') {
    c('peak_value', 'trough_value')
  } else {
    c('diff_mesor', 'diff_peak_trough_amp', 'diff_rms_amp')}#, 'rms_diff_rhy')}
  measCols = intersect(measCols, colnames(posteriorStats))

  d1 = data.table::melt(
    posteriorStats, id.vars = idCols, measure.vars = measCols,
    variable.name = varName, variable.factor = FALSE)
  getInterval = if (method == 'eti') getEti else getHdi
  d2 = d1[, getInterval(value, mass), by = byCols]

  # amp gets special treatment
  if (statType == 'rhy') {
    .SD = NULL
    byCols = intersect(c('feature', 'cond'), colnames(posteriorStats))
    ampCols = intersect(c('peak_trough_amp', 'rms_amp'), colnames(posteriorStats))

    d3 = posteriorStats[, lapply(.SD, getSignedAmp, peak_phase, period),
                        by = byCols, .SDcols = ampCols]
    d4 = data.table::melt(
      d3, id.vars = byCols, variable.name = varName, variable.factor = FALSE)
    d5 = d4[, getInterval(value, mass), by = c(byCols, varName)]
    d2 = rbind(d2, d5)}

  setattr(d2, 'statType', statType)
  setattr(d2, 'period', period)
  setattr(d2, 'mass', mass)
  setattr(d2, 'method', method)
  return(d2)}
