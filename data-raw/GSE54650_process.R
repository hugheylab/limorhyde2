library(data.table)
library(qs)

rawDataPath = file.path('.', 'data-raw')
dataPath = file.path('.', 'data')

d = readRDS(file.path(rawDataPath, 'GSE54650_matrix.rds'))
md = fread(file.path(rawDataPath, 'GSE54650_sample_metadata.csv'))

GSE54650_liver_metadata = md[organ == 'liver']
GSE54650_liver_data = d[, mouse_liver_md$sample]

qsave(GSE54650_liver_metadata, file.path(dataPath, 'GSE54650_liver_metadata.qs'))
qsave(GSE54650_liver_data, file.path(dataPath, 'GSE54650_liver_data.qs'))

fit = getModelFit(y = d, metadata = md, timeColname = 'time')
fit = getPosteriorFit(fit)

rhyStats = getRhythmStats(fit)
fittedVals = getExpectedMeas(fit, times = seq(0, 24, by = 0.5))

#### posterior sampling

fitPs = getPosteriorSamples(fit, nPosteriorSamples = 100)

rhyStatsPs = getRhythmStats(fitPs, fitType = 'posterior_samples')
fittedValsPs = getExpectedMeas(fitPs, times = seq(0, 24, by = 0.5))

qsave(rhyStats, file = 'data/GSE54650_rhystats.qs')
qsave(fittedVals, file = 'data/GSE54650_fittedVals.qs')
qsave(rhyStatsPs, file = 'data/GSE54650_rhystatsPs.qs')
qsave(fittedValsPs, file = 'data/GSE54650_fittedValsPs.qs')

