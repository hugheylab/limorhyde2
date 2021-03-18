

test_that('getModelFit cosinor is functional', {
  library(data.table)

  period = 24
  nKnots = 2

  # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'

  # md0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_md.RDS')
  # d0 = readRDS('/Users/doraobodo/Documents/limorhyde2/tests/testthat/test_zhang_hypoth_data.RDS')


  md0 = readRDS('test_zhang_hypoth_md.RDS')
  d0 = readRDS('test_zhang_hypoth_data.RDS')

  goi = c('22139')#, '69642')

  fit = getModelFit(d0, md0, period,nKnots,
                     timeColname = 'time', conditionsColname = NULL)

  # m = getRhythmAsh(fit)

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



