times = rep(seq(0, 20, 4), 2)
timeColname = 'zt'


test_that('getShifts', {
  period = 18

  s = getShifts(nShifts = 2L, nKnots = 2L, degree = 0L, period = period)
  expect_equal(s, 0)

  s = getShifts(nShifts = 3L, nKnots = 3L, degree = 3L, period = period)
  expect_equal(s, c(0, 1.5, 3))
})


test_that('getMetadata', {
  m1 = data.table(zt = times)
  m1[, con := rep(c('a', 'b'), each = .N / 2)]
  m1[, b1 := rep_len(1:3, .N)]
  m1[, b2 := pi]
  m1[, sample_id := paste0('sample_', seq_len(.N))]

  m = getMetadata(m1, timeColname, condColname = NULL, covarColnames = NULL)
  expect_equal(m, m1[, .(time = zt)])

  m = getMetadata(m1, timeColname, condColname = 'con', covarColnames = NULL)
  expect_equal(m, m1[, .(time = zt, cond = factor(con))])

  m = getMetadata(
    m1, timeColname, condColname = NULL, covarColnames = c('b1', 'b2'))
  expect_equal(m, m1[, .(time = zt, covar_b1 = b1, covar_b2 = b2)])
})


test_that('getBasis', {
  period = 24

  b = getBasis(
    times, period = period, nKnots = 2L, degree = 0L, intercept = TRUE)
  expect_equal(nrow(b), length(times))
  expect_equal(ncol(b), 3L)
  expect_equal(colnames(b), c('intercept', 'basis1', 'basis2'))
  expect_equal(b[1:3, 'basis1'], c(1, 0.5, -0.5))

  nKnots = 4L
  b = getBasis(
    times, period = period, nKnots = nKnots, degree = 3L, intercept = FALSE)
  expect_equal(ncol(b), nKnots)
  expect_equal(colnames(b), paste0('basis', 1:nKnots))
  expect_equal(b[7:8, 'basis1'], c(-0.2, -0.2))
})


test_that('getDesign', {
  m1 = data.table(zt = seq(0, 2.5, 0.5),
                  con = rep(c('a', 'b', 'c'), each = 2L),
                  bat = rep(c('swan', 'elephant'), 3L))

  m = getMetadata(m1, timeColname, condColname = 'con', covarColnames = 'bat')
  d = getDesign(m, period = 1, nKnots = 2L, degree = 0L)

  expect_equal(dim(d), c(nrow(m1), 10L))
  expect_equal(colnames(d)[7:8], c('condb:basis2', 'condc:basis1'))
  expect_equal(unname(d[, 6L]), c(0, 0, 1, -1, 0, 0))
})


test_that('getNumKnotCondCovar', {
  nKnots = 3L
  cols = c('(Intercept)', 'condko', paste0('condwt:basis', 1:nKnots),
           paste0('condko:basis', 1:nKnots), 'covar_batch')
  n = getNumKnotCondCovar(cols)
  expect_equal(n, 3:1)
})
