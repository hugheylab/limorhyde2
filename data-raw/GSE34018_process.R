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
qsave(rhyStats, file = 'data/GSE34018_rhystats.qs')

#### posterior sampling

fitPs = getPosteriorSamples(fit, nPosteriorSamples = 200)
qsave(fitPs, file = 'data/GSE34018_fitPs.qs')