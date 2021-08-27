# x = readRDS('tests/testthat/GSE34018_expression_data.rds')
# metadata = data.table::fread('tests/testthat/GSE34018_sample_metadata.csv')
#
# m = rbind(metadata, metadata[seq(1, .N, 2)])
# m[(.N / 3 * 2 + 1):.N, cond := 'synthetic']
# y = x[seq(1, nrow(x), 200), m$sample]
#
# qs::qsave(y, 'tests/testthat/GSE34018_test_data.qs')
# data.table::fwrite(m, 'tests/testthat/GSE34018_test_metadata.csv')

########################################

foreach::registerDoSEQ()
m = data.table::fread('GSE34018_test_metadata.csv')

timeColname = 'time_test'
data.table::setnames(m, 'time', timeColname)

condColname = 'cond'
conds = data.table(lab = c('wild-type', 'knockout', 'synthetic'),
                   lev = c('wt', 'ko', 'syn'))
set(m, j = condColname, value = factor(m[[condColname]], conds$lab, conds$lev))

m[, batch := rep(c('a', 'b'), length.out = .N)]

y = qs::qread('GSE34018_test_data.qs')
y = y[, m$sample]


test_that('getModelFit', {
  period = 24
  nKnots = 3L

  id = 1
  fitObs = getModelFit(y, m, period, nKnots, timeColname, keepLmFits = TRUE)
  # qs::qsave(fitObs, sprintf('model_fit_%d.qs', id))
  fitExp = qs::qread(sprintf('model_fit_%d.qs', id))
  expect_equal(fitObs, fitExp)

  id = 2
  fitObs = getModelFit(y, m, period, nKnots, timeColname, condColname)
  # qs::qsave(fitObs, sprintf('model_fit_%d.qs', id))
  fitExp = qs::qread(sprintf('model_fit_%d.qs', id))
  expect_equal(fitObs, fitExp)

  id = 3
  fitObs = getModelFit(y, m, period, nKnots, timeColname, condColname,
                       covarColnames = 'batch')
  # qs::qsave(fitObs, sprintf('model_fit_%d.qs', id))
  fitExp = qs::qread(sprintf('model_fit_%d.qs', id))
  expect_equal(fitObs, fitExp)

  id = 4
  fitObs = getModelFit(y, m, period, nKnots = NULL, timeColname = timeColname)
  # qs::qsave(fitObs, sprintf('model_fit_%d.qs', id))
  fitExp = qs::qread(sprintf('model_fit_%d.qs', id))
  expect_equal(fitObs, fitExp)

  expect_error(getModelFit(y, m[-1L], period, nKnots, timeColname))
  expect_error(getModelFit(y, m, -1, nKnots, timeColname))
  expect_error(getModelFit(y, m, period, 2L, timeColname))
  expect_error(getModelFit(y, m, period, nKnots))
})


test_that('getPosteriorFit', {
  id = 1
  fit = qs::qread(sprintf('model_fit_%d.qs', id))
  fitObs = getPosteriorFit(fit)
  # qs::qsave(fitObs, sprintf('posterior_fit_%d.qs', id))
  fitExp = qs::qread(sprintf('posterior_fit_%d.qs', id))

  expect_equal(fitObs, fitExp)
  expect_error(getPosteriorFit(fitObs))
  # expect_equal(getPosteriorFit(fitObs, overwrite = TRUE), fitExp)

  # id = 1
  # fit = qs::qread(sprintf('model_fit_%d.qs', id))
  # fitObs = getPosteriorFit(fit, covMethod = 'canonical')
  # # qs::qsave(fitObs, sprintf('posterior_fit_canon_%d.qs', id))
  # fitExp = qs::qread(sprintf('posterior_fit_canon_%d.qs', id))
  #
  # expect_equal(fitObs, fitExp)
})


test_that('getRhythmStats', {
  id = 1
  fit = qs::qread(sprintf('model_fit_%d.qs', id))
  statsObs = getRhythmStats(fit, 'raw')
  # qs::qsave(statsObs, sprintf('rhy_stats_raw_%d.qs', id))
  statsExp = qs::qread(sprintf('rhy_stats_raw_%d.qs', id))

  expect_equal(statsObs, statsExp)
  expect_error(getRhythmStats(fit))

  id = 2
  fit = qs::qread(sprintf('posterior_fit_%d.qs', id))
  statsObs = getRhythmStats(fit)
  # qs::qsave(statsObs, sprintf('rhy_stats_post_%d.qs', id))
  statsExp = qs::qread(sprintf('rhy_stats_post_%d.qs', id))

  expect_equal(statsObs, statsExp)
  expect_equal(getRhythmStats(fit, features = '11287'),
               statsExp[feature == '11287'])

  fit = qs::qread(sprintf('posterior_samples_%d.qs', id))
  statsObs = getRhythmStats(fit, 'posterior_samples', features = 1:2)
  # qs::qsave(statsObs, sprintf('rhy_stats_samps_%d.qs', id))
  statsExp = qs::qread(sprintf('rhy_stats_samps_%d.qs', id))

  expect_equal(statsObs, statsExp)
})


test_that('getDiffRhythmStats', {
  id = 2
  fit = qs::qread(sprintf('posterior_fit_%d.qs', id))
  rhyStats = qs::qread(sprintf('rhy_stats_post_%d.qs', id))
  statsObs = getDiffRhythmStats(fit, rhyStats)
  # qs::qsave(statsObs, sprintf('diff_rhy_stats_post_%d.qs', id))
  statsExp = qs::qread(sprintf('diff_rhy_stats_post_%d.qs', id))

  expect_equal(statsObs, statsExp)
  expect_equal(levels(statsObs$cond1), fit$condLevels)
  expect_equal(levels(statsObs$cond2), fit$condLevels)

  fit = qs::qread(sprintf('posterior_samples_%d.qs', id))
  rhyStats = qs::qread(sprintf('rhy_stats_samps_%d.qs', id))
  statsObs = getDiffRhythmStats(fit, rhyStats, conds$lev[1:2])
  # qs::qsave(statsObs, sprintf('diff_rhy_stats_samps12_%d.qs', id))
  statsExp = qs::qread(sprintf('diff_rhy_stats_samps12_%d.qs', id))

  expect_equal(statsObs, statsExp)

  statsObs = getDiffRhythmStats(fit, rhyStats, conds$lev[3:2])
  # qs::qsave(statsObs, sprintf('diff_rhy_stats_samps32_%d.qs', id))
  statsExp = qs::qread(sprintf('diff_rhy_stats_samps32_%d.qs', id))

  expect_equal(statsObs, statsExp)

  id = 1
  fit = qs::qread(sprintf('model_fit_%d.qs', id))
  rhyStats = qs::qread(sprintf('rhy_stats_raw_%d.qs', id))

  expect_error(getDiffRhythmStats(fit, rhyStats, conds$lev[1:2]))
})


test_that('getExpectedMeas', {
  times = 22:27
  features = 1:5

  id = 2
  fit = qs::qread(sprintf('posterior_samples_%d.qs', id))
  measObs = getExpectedMeas(
    fit, times = times, fitType = 'posterior_samples', features = features)
  # qs::qsave(measObs, sprintf('expected_meas_%d.qs', id))
  measExp = qs::qread(sprintf('expected_meas_%d.qs', id))

  expect_equal(measObs, measExp)

  id = 3
  fit = qs::qread(sprintf('model_fit_%d.qs', id))
  measObs = getExpectedMeas(
    fit, times = times, fitType = 'raw', features = features)
  # qs::qsave(measObs, sprintf('expected_meas_%d.qs', id))
  measExp = qs::qread(sprintf('expected_meas_%d.qs', id))

  expect_equal(measObs, measExp)
})


test_that('getPosteriorSamples', {
  id = 2
  fitObs = qs::qread(sprintf('posterior_fit_%d.qs', id))
  fitObs = getPosteriorSamples(fitObs, nPosteriorSamples = 10)
  # qs::qsave(fitObs, sprintf('posterior_samples_%d.qs', id))
  fitExp = qs::qread(sprintf('posterior_samples_%d.qs', id))

  expect_equal(fitObs, fitExp)
  expect_error(getPosteriorSamples(fitObs))
  # expect_equal(getPosteriorSamples(
  #   fitObs, nPosteriorSamples = 10, overwrite = TRUE), fitExp)
})


test_that('getExpectedMeasIntervals', {
  id = 2
  meas = qs::qread(sprintf('expected_meas_%d.qs', id))
  intsObs = getExpectedMeasIntervals(meas)
  # qs::qsave(intsObs, sprintf('meas_ints_%d.qs', id))
  intsExp = qs::qread(sprintf('meas_ints_%d.qs', id))

  expect_equal(intsObs, intsExp)

  id = 3
  meas = qs::qread(sprintf('expected_meas_%d.qs', id))
  expect_error(getExpectedMeasIntervals(meas))
})


test_that('getStatsIntervals', {
  id = 2
  rhyStats = qs::qread(sprintf('rhy_stats_samps_%d.qs', id))
  intsObs = getStatsIntervals(rhyStats)
  # qs::qsave(intsObs, sprintf('stats_ints_rhy_%d.qs', id))
  intsExp = qs::qread(sprintf('stats_ints_rhy_%d.qs', id))

  expect_equal(intsObs, intsExp)

  diffRhyStats = qs::qread(sprintf('diff_rhy_stats_samps12_%d.qs', id))
  intsObs = getStatsIntervals(diffRhyStats)
  # qs::qsave(intsObs, sprintf('stats_ints_diff_%d.qs', id))
  intsExp = qs::qread(sprintf('stats_ints_diff_%d.qs', id))

  expect_equal(intsObs, intsExp)
  expect_error(getStatsIntervals(diffRhyStats, mass = 90))
})


test_that('mergeMeasMeta', {
  nSamps = 2L
  sampleColname = 'sample_z'

  metadata = data.table(sample = paste0('sample_', 1:nSamps),
                        time = seq(0, 20, length.out = nSamps))
  data.table::setnames(metadata, 'sample', sampleColname)

  y = matrix(1:8, ncol = nSamps)
  colnames(y) = metadata[[sampleColname]]
  rownames(y) = paste0('feature_', 1:nrow(y))

  dObs = mergeMeasMeta(
    y, metadata, features = rownames(y)[1L], sampleColname = sampleColname)

  dExp = data.table(
    metadata, feature = rownames(y)[1L], meas = y[1L, ], key = sampleColname)

  expect_equal(dObs, dExp)
})
