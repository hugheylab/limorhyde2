library(data.table)
library(qs)
library(limorhyde2)

rawDataPath = file.path('.', 'data-raw')
dataPath = file.path('.', 'data')

d = readRDS(file.path(rawDataPath, 'GSE34018_expression_data.rds'))
md = fread(file.path(rawDataPath, 'GSE34018_sample_metadata.csv'))

qsave(md, file.path(dataPath, 'GSE34018_metadata.qs'))
qsave(d, file.path(dataPath, 'GSE34018_data.qs'))

fit = getModelFit(y = d, metadata = md, timeColname = 'time',
                  condColname = 'cond')
fit = getPosteriorFit(fit)

rhyStats = getRhythmStats(fit)
fittedVals = getExpectedMeas(fit, times = seq(0, 24, by = 0.5))

#### posterior sampling

fitPs = getPosteriorSamples(fit, nPosteriorSamples = 100)

rhyStatsPs = getRhythmStats(fitPs, fitType = 'posterior_samples')
fittedValsPs = getExpectedMeas(fitPs, times = seq(0, 24, by = 0.5))

qsave(rhyStats, file = 'data/GSE34018_rhystats.qs')
qsave(fittedVals, file = 'data/GSE34018_fittedVals.qs')
qsave(rhyStatsPs, file = 'data/GSE34018_rhystatsPs.qs')
qsave(fittedValsPs, file = 'data/GSE34018_fittedValsPs.qs')
