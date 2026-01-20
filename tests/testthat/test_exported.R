y = GSE34018$y
m = GSE34018$metadata

m = rbind(m, m[seq(1, .N, 2)])
m[(.N / 3 * 2 + 1):.N, cond := 'synthetic']
y = y[, m$sample]
m[, sample := paste0(sample, '_', seq_len(.N))]
colnames(y) = m$sample

timeColname = 'time_test'
data.table::setnames(m, 'time', timeColname)

condColname = 'cond'
conds = data.table(lab = c('wild-type', 'knockout', 'synthetic'),
                   lev = c('wt', 'ko', 'syn'))
set(m, j = condColname, value = factor(m[[condColname]], conds$lab, conds$lev))

m[, batch := rep(c('a', 'b'), length.out = .N)]


test_that('getModelFit', {
  period = 24
  nKnots = 3L

  id = 1
  fitObs = getModelFit(
    y, m, period, nKnots, timeColname = timeColname, keepLmFits = TRUE)
  fitExp = snapshot(fitObs, file.path(dataDir, glue('model_fit_{id}.rds')))
  expect_equal(fitObs, fitExp)

  id = 2
  fitObs = getModelFit(
    y, m, period, nKnots, timeColname = timeColname, condColname = condColname)
  fitExp = snapshot(fitObs, file.path(dataDir, glue('model_fit_{id}.rds')))
  expect_equal(fitObs, fitExp)

  id = 3
  fitObs = getModelFit(
    y, m, period, nKnots, timeColname = timeColname, condColname = condColname,
    covarColnames = 'batch')
  fitExp = snapshot(fitObs, file.path(dataDir, glue('model_fit_{id}.rds')))
  expect_equal(fitObs, fitExp)

  id = 4
  fitObs = getModelFit(y, m, period, sinusoid = TRUE, timeColname = timeColname)
  fitExp = snapshot(fitObs, file.path(dataDir, glue('model_fit_{id}.rds')))
  expect_equal(fitObs, fitExp)

  expect_error(getModelFit(y, m[-1L], period, nKnots, timeColname = timeColname))
  expect_error(getModelFit(y, m, -1, nKnots, timeColname = timeColname))
  expect_error(getModelFit(y, m, period, nKnots))
})


test_that('getPosteriorFit', {
  id = 2
  fit = readRDS(file.path(dataDir, glue('model_fit_{id}.rds')))
  fitObs = getPosteriorFit(fit)
  fitExp = snapshot(
    fitObs, file.path(dataDir, glue('posterior_fit_{id}.rds')))

  expect_equal(fitObs, fitExp, tolerance = 1e-4) # avoid intermittent failures
  expect_error(getPosteriorFit(fitObs))
  # expect_equal(getPosteriorFit(fitObs, overwrite = TRUE), fitExp)

  # id = 1
  # fit = readRDS(glue('model_fit_{id}.rds'))
  # fitObs = getPosteriorFit(fit, covMethod = 'canonical')
  # # qs::qsave(fitObs, glue('posterior_fit_canon_{id}.rds'))
  # fitExp = readRDS(glue('posterior_fit_canon_{id}.rds'))
  #
  # expect_equal(fitObs, fitExp)
})


test_that('getPosteriorSamples', {
  id = 2
  fitObs = readRDS(file.path(dataDir, glue('posterior_fit_{id}.rds')))
  fitObs = getPosteriorSamples(fitObs, nPosteriorSamples = 10)
  fitExp = snapshot(
    fitObs, file.path(dataDir, glue('posterior_samples_{id}.rds')))

  expect_equal(fitObs, fitExp)
  expect_error(getPosteriorSamples(fitObs))
  # expect_equal(getPosteriorSamples(
  #   fitObs, nPosteriorSamples = 10, overwrite = TRUE), fitExp)
})


test_that('getRhythmStats', {
  id = 1
  fit = readRDS(file.path(dataDir, glue('model_fit_{id}.rds')))
  statsObs = getRhythmStats(fit, 'raw')
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('rhy_stats_raw_{id}.rds')))

  expect_equal(statsObs, statsExp)
  expect_error(getRhythmStats(fit))

  id = 2
  fit = readRDS(file.path(dataDir, glue('posterior_fit_{id}.rds')))
  statsObs = getRhythmStats(fit, rms = TRUE)
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('rhy_stats_post_{id}.rds')))

  expect_equal(statsObs, statsExp)
  expect_equal(getRhythmStats(fit, features = '11287', rms = TRUE),
               statsExp[feature == '11287'])

  fit = readRDS(file.path(dataDir, glue('posterior_samples_{id}.rds')))
  statsObs = getRhythmStats(fit, 'posterior_samples', features = 1:2)
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('rhy_stats_samps_{id}.rds')))

  expect_equal(statsObs, statsExp)
})


test_that('getDiffRhythmStats', {
  id = 2
  fit = readRDS(file.path(dataDir, glue('posterior_fit_{id}.rds')))
  rhyStats = readRDS(file.path(dataDir, glue('rhy_stats_post_{id}.rds')))
  statsObs = getDiffRhythmStats(fit, rhyStats)
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('diff_rhy_stats_post_{id}.rds')))

  expect_equal(statsObs, statsExp)

  fit = readRDS(file.path(dataDir, glue('posterior_samples_{id}.rds')))
  rhyStats = readRDS(file.path(dataDir, glue('rhy_stats_samps_{id}.rds')))
  statsObs = getDiffRhythmStats(fit, rhyStats, conds$lev[1:2])
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('diff_rhy_stats_samps12_{id}.rds')))

  expect_equal(statsObs, statsExp)

  statsObs = getDiffRhythmStats(fit, rhyStats, conds$lev[3:2])
  statsExp = snapshot(
    statsObs, file.path(dataDir, glue('diff_rhy_stats_samps32_{id}.rds')))

  expect_equal(statsObs, statsExp)

  id = 1
  fit = readRDS(file.path(dataDir, glue('model_fit_{id}.rds')))
  rhyStats = readRDS(file.path(dataDir, glue('rhy_stats_raw_{id}.rds')))

  expect_error(getDiffRhythmStats(fit, rhyStats, conds$lev[1:2]))
})


test_that('getExpectedMeas', {
  times = 22:27
  features = 1:5

  id = 2
  fit = readRDS(file.path(dataDir, glue('posterior_samples_{id}.rds')))
  measObs = getExpectedMeas(
    fit, times = times, fitType = 'posterior_samples', features = features)
  measExp = snapshot(
    measObs, file.path(dataDir, glue('expected_meas_{id}.rds')))

  expect_equal(measObs, measExp)

  id = 3
  fit = readRDS(file.path(dataDir, glue('model_fit_{id}.rds')))
  measObs = getExpectedMeas(
    fit, times = times, fitType = 'raw', features = features)
  measExp = snapshot(
    measObs, file.path(dataDir, glue('expected_meas_{id}.rds')))

  expect_equal(measObs, measExp)
})


test_that('getExpectedMeasIntervals', {
  id = 2
  meas = readRDS(file.path(dataDir, glue('expected_meas_{id}.rds')))
  intsObs = getExpectedMeasIntervals(meas)
  intsExp = snapshot(
    intsObs, file.path(dataDir, glue('meas_ints_{id}.rds')))

  expect_equal(intsObs, intsExp)

  id = 3
  meas = readRDS(file.path(dataDir, glue('expected_meas_{id}.rds')))
  expect_error(getExpectedMeasIntervals(meas))
})


test_that('getStatsIntervals', {
  id = 2
  rhyStats = readRDS(file.path(dataDir, glue('rhy_stats_samps_{id}.rds')))
  intsObs = getStatsIntervals(rhyStats)
  intsExp = snapshot(
    intsObs, file.path(dataDir, glue('stats_ints_rhy_{id}.rds')))

  expect_equal(intsObs, intsExp)

  diffRhyStats = readRDS(
    file.path(dataDir, glue('diff_rhy_stats_samps12_{id}.rds')))
  intsObs = getStatsIntervals(diffRhyStats)
  intsExp = snapshot(
    intsObs, file.path(dataDir, glue('stats_ints_diff_{id}.rds')))

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
  rownames(y) = paste0('feature_', seq_len(nrow(y)))

  dObs = mergeMeasMeta(
    y, metadata, features = rownames(y)[1L], sampleColname = sampleColname)

  dExp = data.table(
    metadata, feature = rownames(y)[1L], meas = y[1L, ], key = sampleColname)

  expect_equal(dObs, dExp)
})
