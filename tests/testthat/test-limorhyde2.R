context('Exported functions')

test_that('getModelFit cosinor is functional', {
  library(data.table)

  period = 24
  nKnots = 4

  # path = '/Users/doraobodo/Documents/limorhyde2/tests/testthat'

  d0 = read.csv('test_limorhyde2_one_cond_example_data.csv', row.names = 1)
  md0 = read.csv('test_limorhyde2_one_cond_example_md.csv', row.names = 1)


  fit = getModelFit(d0, md0, period,nKnots,
                    timeColname = 'time', conditionsColname = NULL)
  fitC = as.data.table(fit$coefficients, rownames = TRUE)
  # fwrite(fitC, file = 'cosinor_fit_coefs_test_output.csv')

  fitE = fread('test_limorhyde2_one_cond_example_cosinor_fit_coefs.csv')

  fitEqual = all.equal(fitC, fitE, check.attributes = FALSE)
  # write(fitEqual, file = 'all_eq_test1.txt')
  expect_true(fitEqual)

})