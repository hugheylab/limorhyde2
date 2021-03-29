

test_that('getModelFit cosinor is functional', {
  library(data.table)

  period = 24
  nKnots = 4

  # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'

  md0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_md.RDS')
  d0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_data.RDS')


  md0 = readRDS('test_zhang_hypoth_md.RDS')
  d0 = readRDS('test_zhang_hypoth_data.RDS')

  goi = c('22139')#, '69642')

  fit = getModelFit(d0[1:3,], md0, period,nKnots = nKnots,
                     timeColname = 'time', conditionsColname = NULL)


  m = getRhythmAsh(fit, 'canonical', npc = nKnots+1)
  #
  cStats = getRhythmStats(fit$coefficients)
  mStats = getRhythmStats(m)
  #
#
#   fitC = getModelFit(d0, md0, period,nKnots,
#                     timeColname = 'time', conditionsColname = NULL)
#
#   mashC = getRhythmAsh(fitC)
#
#   fitR = getModelFit(d0, md0, period,nKnots,
#                      timeColname = 'time', conditionsColname = NULL)
#
#   mashR = getRhythmAsh(fitR)
#
#   fit$coefficients[goi,]
#   m[goi,]
#
#   fitC$coefficients[goi,]
#   mashC[goi,]
#
#   fitR$coefficients[goi,]
#   mashR[goi,]

})

test_that('getRhythmStats spline is functional for two conditions', {
  library(data.table)

  period = 24
  nKnots = 4

  # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'

  md0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_limorhyde2_two_cond_example_md.rds')
  d0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_limorhyde2_two_cond_example_data.rds')

  fit = getModelFit(d0[1:100,], md0, period,nKnots,
                    timeColname = 'time', conditionsColname = 'cond')

  c = fit$coefficients
  m = getRhythmAsh(fit, 'data-driven', npc = nKnots-1)

  cStats = getRhythmStats(fit$coefficients)
  mStats = getRhythmStats(m)

  cDiffStats = getDiffRhythmStats(cStats, condIds = c("wt",'ko'))
  mDiffStats = getDiffRhythmStats(mStats, condIds = c("wt",'ko'))


  ### check equality of stats for each condition

  # # WT
  # md0 = md0[cond == 'wt']
  # d0 = d0[, md0$sample]
  #
  # fit = getModelFit(d0, md0, period,nKnots,
  #                   timeColname = 'time', conditionsColname = NULL)
  #
  # cWT = fit$coefficients
  # mWT = getRhythmAsh(fit, 'data-driven', npc = nKnots-1)
  #
  # cStatsWT = getRhythmStats(cWT)
  # mStatsWT = getRhythmStats(mWT)
  #
  # # KO
  # md0 = md0[cond == 'ko']
  # d0 = d0[, md0$sample]
  #
  # fit = getModelFit(d0, md0, period,nKnots,
  #                   timeColname = 'time', conditionsColname = NULL)
  #
  # cKO = fit$coefficients
  # mKO = getRhythmAsh(fit, 'data-driven', npc = nKnots-1)
  #
  # cStatsKO = getRhythmStats(cKO)
  # mStatsKO = getRhythmStats(mKO)
  #
  #
  # cDiff = cKO - cWT
  # mDiff = mKO - mWT
  #
  # cStatsSub = as.matrix(cStatsKO, rownames = 'feature') - as.matrix(cStatsWT, rownames = 'feature')
  # mStatsSub = mStatsKO - mStatsWT
#
})


