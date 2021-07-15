library('data.table')
library('limorhyde2')
library('qs')

doParallel::registerDoParallel()
rawDataDir = 'data-raw'
dataDir = file.path('inst', 'extdata')
nGenes = 500L

metadata = fread(file.path(rawDataDir, 'GSE54650_sample_metadata.csv'))
y = qread(file.path(rawDataDir, 'GSE54650_expression_data.qs'))

metadata = metadata[tissue == 'liver', .(sample, time)]

idx = c(which(rownames(y) %in% c('13088', '13170', '13869')),
        round(seq(1, nrow(y), length.out = nGenes - 3)))
y = y[idx, metadata$sample]

qsave(metadata, file.path(dataDir, 'GSE54650_liver_metadata.qs'))
qsave(y, file.path(dataDir, 'GSE54650_liver_data.qs'))

# fit = getModelFit(y, metadata)
# fit = getPosteriorFit(fit)
#
# rhyStats = getRhythmStats(fit)
# qsave(rhyStats, file.path(dataDir, 'GSE54650_liver_rhy_stats.qs'))
