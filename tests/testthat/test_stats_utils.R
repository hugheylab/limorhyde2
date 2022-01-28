test_that('checkFitType', {
  fit = list(mashCoefficients = 1L)
  expect_invisible(checkFitType(fit, 'posterior_mean'))
  expect_error(checkFitType(fit, 'posterior_samples'))

  fit = list(mashPosteriorSamples = 1L)
  expect_error(checkFitType(fit, 'posterior_mean'))
  expect_invisible(checkFitType(fit, 'posterior_samples'))
})


test_that('getCoefArray', {
  n = c(7L, 3L, 5L)
  nNames = list(paste0('feature_', 1:n[1L]), paste0('coef_', 1:n[2L]), NULL)
  d = data.table(fitType = c('raw', 'posterior_mean', 'posterior_samples'),
                 value = 0:2)

  fit = list(
    coefficients = matrix(d$value[1L], nrow = n[1L], ncol = n[2L]),
    mashCoefficients = matrix(d$value[2L], nrow = n[1L], ncol = n[2L]),
    mashPosteriorSamples = array(d$value[3L], dim = n, dimnames = nNames))
  for (i in 1:2) dimnames(fit[[i]]) = nNames[1:2]

  for (i in 1:2) {
    coArray = getCoefArray(fit, d$fitType[i])
    expect_equal(coArray[1L], d$value[i])
    expect_equal(dim(coArray), c(n[1:2], 1L))
    expect_equal(dimnames(coArray), nNames)}

  i = 3
  coArray = getCoefArray(fit, d$fitType[i])
  expect_equal(coArray[1L], d$value[i])
  expect_equal(dim(coArray), n)
  expect_equal(dimnames(coArray), nNames)
})


test_that('getOptima', {
  tr = seq(0, 2*pi, length.out = 40)
  f = sin
  dObs = getOptima(f, tr)
  dExp = data.table(peak_phase = pi / 2, peak_value = 1,
                    trough_phase = 3 * pi / 2, trough_value = -1)
  expect_equal(dObs, dExp)
})


test_that('getCoefMatOneCond', {
  nConds = 2L
  nKnots = 2L
  nShifts = 2L
  cols = c('(Intercept)', 'condb', paste0('conda:basis', 1:nKnots),
           paste0('condb:basis', 1:nKnots))
  cols = paste0(rep(cols, nShifts),
                rep(paste0('_shift', 1:nShifts), each = length(cols)))

  nFeats = 1L
  condIdx = 1L
  coefMat = matrix(rep(1:length(cols), each = nFeats), ncol = length(cols))
  coefObs = getCoefMatOneCond(coefMat, condIdx, nConds, nKnots, nShifts)
  coefExp = matrix(c(1, 3, 4, 7, 9, 10), nrow = nFeats)
  expect_equal(coefObs, coefExp)

  condIdx = 2L
  coefMat = matrix(rep(1:length(cols), each = nFeats), ncol = length(cols))
  coefObs = getCoefMatOneCond(coefMat, condIdx, nConds, nKnots, nShifts)
  coefExp = matrix(c(3, 5, 6, 15, 11, 12), nrow = nFeats)
  expect_equal(coefObs, coefExp)

  nFeats = 8L
  condIdx = 1L
  coefMat = matrix(rep(1:length(cols), each = nFeats), ncol = length(cols))
  coefObs = getCoefMatOneCond(coefMat, condIdx, nConds, nKnots, nShifts)
  coefExp = matrix(rep(c(1, 3, 4, 7, 9, 10), each = nFeats), nrow = nFeats)
  expect_equal(coefObs, coefExp)

  condIdx = 2L
  coefMat = matrix(rep(1:length(cols), each = nFeats), ncol = length(cols))
  coefObs = getCoefMatOneCond(coefMat, condIdx, nConds, nKnots, nShifts)
  coefExp = matrix(rep(c(3, 5, 6, 15, 11, 12), each = nFeats), nrow = nFeats)
  expect_equal(coefObs, coefExp)
})


test_that('getRmsAmp', {
  f = function(x) sin(x) + 42
  co = c(42, -1.492, 1984)
  period = 2 * pi
  r = getRmsAmp(f, co, period)
  expect_equal(r, sqrt(2) / 2)
})


test_that('centerCircDiff', {
  x = c(-7, -5, -3, 0, 3, 5, 7)
  period = 10
  r = centerCircDiff(x, period)
  expect_equal(r, c(3, 5, -3, 0, 3, 5, -3))
})


test_that('getEti', {
  v = seq(2, 3, 0.01)
  mass = 0.77
  dObs = getEti(v, mass)
  dExp = data.table(lower = v[1L] + (1 - mass) / 2,
                    upper = v[1L] + (1 + mass) / 2)
  expect_equal(dObs, dExp)
})


test_that('getHdi', {
  mass = 0.8
  v = c(rep(0:2, 10), rep(3:5, 80), rep(6:8, 10))
  dObs = getHdi(v, mass)
  dExp = data.table(lower = 2, upper = 5)
  expect_equal(dObs, dExp)
})
