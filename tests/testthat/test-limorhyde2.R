# test_that('getModelFit cosinor is functional', {
#   # library(data.table)
#
#   period = 24
#   nKnots = 4
#
#   # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'
#
#   # md0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_md.RDS')
#   # d0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_data.RDS')
#
#
#   md0 = readRDS('test_zhang_hypoth_md.RDS')
#   d0 = readRDS('test_zhang_hypoth_data.RDS')
#
#   goi = c('22139')#, '69642')
#
#   fit = getModelFit(d0[1:1000, ], md0, period, nKnots = nKnots,
#                      timeColname = 'time', condColname = NULL)
#
#   # fit2 = getExpectedMeas(fit, times = 3, fitType = 'raw')
#
#   m = getPosteriorFit(fit, 'canonical', npc = nKnots + 1)
#   #
#   # cStats = getRhythmStats(fit$coefficients)
#   mStats = getRhythmStats(m, features = 1:10)
#   #
# #
# #   fitC = getModelFit(d0, md0, period,nKnots,
# #                     timeColname = 'time', conditionsColname = NULL)
# #
# #   mashC = getPosteriorFit(fitC)
# #
# #   fitR = getModelFit(d0, md0, period,nKnots,
# #                      timeColname = 'time', conditionsColname = NULL)
# #
# #   mashR = getPosteriorFit(fitR)
# #
# #   fit$coefficients[goi,]
# #   m[goi,]
# #
# #   fitC$coefficients[goi,]
# #   mashC[goi,]
# #
# #   fitR$coefficients[goi,]
# #   mashR[goi,]
#
# })
#
# test_that('getRhythmStats spline is functional for two conditions', {
#   # library(data.table)
#
#   period = 24
#   nKnots = 4
#
#   # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'
#
#   # md0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_limorhyde2_two_cond_example_md.rds')
#   # d0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_limorhyde2_two_cond_example_data.rds')
#
#   fit = getModelFit(d0[1:1000, ], md0, period, nKnots, timeColname = 'time',
#                     condColname = 'cond')
#
#   # c = fit$coefficients
#   m = getPosteriorFit(fit, 'data-driven', npc = nKnots - 1)
#
#   # cStats = getRhythmStats(fit$coefficients)
#   mStats = getRhythmStats(m, features = 1:10)
#
#   # cDiffStats = getDiffRhythmStats(cStats, condLevels = c('wt', 'ko'))
#   mDiffStats = getDiffRhythmStats(mStats, condLevels = c('wt', 'ko'))
#
#
#   ### check equality of stats for each condition
#
#   # # WT
#   # md0 = md0[cond == 'wt']
#   # d0 = d0[, md0$sample]
#   #
#   # fit = getModelFit(d0, md0, period,nKnots,
#   #                   timeColname = 'time', conditionsColname = NULL)
#   #
#   # cWT = fit$coefficients
#   # mWT = getPosteriorFit(fit, 'data-driven', npc = nKnots-1)
#   #
#   # cStatsWT = getRhythmStats(cWT)
#   # mStatsWT = getRhythmStats(mWT)
#   #
#   # # KO
#   # md0 = md0[cond == 'ko']
#   # d0 = d0[, md0$sample]
#   #
#   # fit = getModelFit(d0, md0, period,nKnots,
#   #                   timeColname = 'time', conditionsColname = NULL)
#   #
#   # cKO = fit$coefficients
#   # mKO = getPosteriorFit(fit, 'data-driven', npc = nKnots-1)
#   #
#   # cStatsKO = getRhythmStats(cKO)
#   # mStatsKO = getRhythmStats(mKO)
#   #
#   #
#   # cDiff = cKO - cWT
#   # mDiff = mKO - mWT
#   #
#   # cStatsSub = as.matrix(cStatsKO, rownames = 'feature') - as.matrix(cStatsWT, rownames = 'feature')
#   # mStatsSub = mStatsKO - mStatsWT
# #
# })
#
#
