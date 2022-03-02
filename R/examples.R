# https://stackoverflow.com/questions/59327170/can-you-share-examples-between-functions-with-roxygen2

examples1 = function() {
  ex = "
@examples
library('data.table')

# rhythmicity in one condition
y = GSE54650$y
metadata = GSE54650$metadata

fit = getModelFit(y, metadata)
fit = getPosteriorFit(fit)
rhyStats = getRhythmStats(fit, features = c('13170', '13869'))

# rhythmicity and differential rhythmicity in multiple conditions
y = GSE34018$y
metadata = GSE34018$metadata

fit = getModelFit(y, metadata, nKnots = 3L, condColname = 'cond')
fit = getPosteriorFit(fit)
rhyStats = getRhythmStats(fit, features = c('13170', '12686'))
diffRhyStats = getDiffRhythmStats(fit, rhyStats)
"
  return(strsplit(ex, split = '\n')[[1L]])}


examples2 = function() {
  ex = "
@examples
library('data.table')

y = GSE54650$y
metadata = GSE54650$metadata

fit = getModelFit(y, metadata)
fit = getPosteriorFit(fit)
fit = getPosteriorSamples(fit, nPosteriorSamples = 10L)

rhyStatsSamps = getRhythmStats(
  fit, features = c('13170', '13869'), fitType = 'posterior_samples')
rhyStatsInts = getStatsIntervals(rhyStatsSamps)
"
  return(strsplit(ex, split = '\n')[[1L]])}


examples3 = function() {
  ex = "
@examples
library('data.table')

y = GSE34018$y
metadata = GSE34018$metadata

fit = getModelFit(y, metadata)
fit = getPosteriorFit(fit)

measObs = mergeMeasMeta(y, metadata, features = c('13170', '12686'))
measFitMean = getExpectedMeas(
  fit, times = seq(0, 24, 0.5), features = c('13170', '12686'))
"
  return(strsplit(ex, split = '\n')[[1L]])}


examples4 = function() {
  ex = "
@examples
library('data.table')

y = GSE34018$y
metadata = GSE34018$metadata

fit = getModelFit(y, metadata)
fit = getPosteriorFit(fit)
fit = getPosteriorSamples(fit, nPosteriorSamples = 10L)

measFitSamps = getExpectedMeas(
  fit, times = seq(0, 24, 0.5), fitType = 'posterior_samples',
  features = c('13170', '12686'))
measFitInts = getExpectedMeasIntervals(measFitSamps)
"
  return(strsplit(ex, split = '\n')[[1L]])}

