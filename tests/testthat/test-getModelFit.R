# test_that('getBasis intercept argument is binary and returns column in a matrix ', {
#   x = seq(0, 24, 2)
#   period = 12
#   nKnots = 3
#   b1 = getBasis(x, period, nKnots = 2, intercept = TRUE)
#   b0 = getBasis(x, period, nKnots=nKnots, intercept = FALSE)
#
#   expect_error(getBasis(x, intercept = 0))
#   expect_type(b1, 'double')
#   expect_identical(colnames(b1)[1], 'intercept')
#   })
#
# test_that('getBasis period and nKnots values return the correct number of rows and columns', {
#
#   x = seq(0, 24, 2)
#   p = 12
#   n = 3
#
#   b = getBasis(x, period = p, nKnots = n, intercept = FALSE)
#
#
#   expect_error(getBasis(x, 'cat'))
#   expect_error(getBasis(x, c('cat', 24)))
#   expect_error(getBasis(x, -2))
#
#   expect_error(getBasis(x, nKnots = 'cat'))
#   expect_error(getBasis(x, nKnots = c('cat', 24)))
#   expect_error(getBasis(x, nKnots = 1))
#
#   expect_equal(ncol(b), n)
#   expect_equal(nrow(b), length(x))
#
#
# })
#
# test_that('getMetadata returns a data.table with specific columns', {
#
#   # 1. time column name has to be specified and hcolumn has to be numerics only
#   period = 12
#   samples = rep(seq(0,period,2.3),each =2)
#   mdLetter = data.table(time = letters[1:length(samples)], sample_id = paste0('sample_', 1:length(samples)))
#   mdNum = data.table(time = samples, sample_id = paste0('sample_', 1:length(samples)),
#                      conds = rep_len(c('wt', 'ko'), length.out = length(samples) ),
#                      misc = c('test'))
#   d = getMetadata(mdNum, 'time', 'conds')
#
#
#   expect_error(getMetadata(mdLetter, 'timed', NULL))
#   expect_error(getMetadata(mdLetter, 1, NULL))
#   expect_error(getMetadata(mdLetter, c('timed', 'd'), NULL))
#   expect_error(getMetadata(mdLetter, 'time', NULL))
#   expect_error(getMetadata(mdNum, 'time', 'conditions'))
#
#   expect_s3_class(getMetadata(mdNum, 'time', NULL), 'data.table')
#   expect_equal(colnames(d), c('time', 'cond'))
#
#
# })
#
# test_that('getModelFit', {
#
#   y = readRDS('simulate_3cond_gene_expression.RDS')
#   #check that sample rownames match sample column names
#   md = data.table(time = rep(seq(0,23,4),3), sample_id = colnames(y), cond = rep(c('wt','ko', 'vv'),
#                   each = ncol(y)/3), misc = paste('Unn_info'))
#
#
#
#
# })
