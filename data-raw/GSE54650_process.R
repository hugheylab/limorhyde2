library(data.table)
library(qs)
library(limorhyde2)

rawDataPath = file.path('.', 'data-raw')
dataPath = file.path('.', 'data')

d = readRDS(file.path(rawDataPath, 'GSE54650_matrix.rds'))
md = fread(file.path(rawDataPath, 'GSE54650_sample_metadata.csv'))

GSE54650_liver_metadata = md[organ == 'liver']
GSE54650_liver_data = d[, GSE54650_liver_metadata$sample]

qsave(GSE54650_liver_metadata, file.path(dataPath, 'GSE54650_liver_metadata.qs'))
qsave(GSE54650_liver_data, file.path(dataPath, 'GSE54650_liver_data.qs'))

fit = getModelFit(y = GSE54650_liver_data, metadata = GSE54650_liver_metadata,
                  timeColname = 'time')
fit = getPosteriorFit(fit)

rhyStats = getRhythmStats(fit) #rhythmicity vignette
qsave(rhyStats, file = 'data/GSE54650_rhystats.qs')

#### posterior sampling

fitPs = getPosteriorSamples(fit, nPosteriorSamples = 200) #diff rhythmicity vignette
qsave(fitPs, file = 'data/GSE54650_fitPs.qs')
