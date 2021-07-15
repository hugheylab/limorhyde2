library('data.table')
library('limorhyde2')
library('qs')

doParallel::registerDoParallel()
rawDataDir = 'data-raw'
dataDir = file.path('inst', 'extdata')
nGenes = 1000L

metadata = fread(file.path(rawDataDir, 'GSE34018_sample_metadata.csv'))
y = qread(file.path(rawDataDir, 'GSE34018_expression_data.qs'))

metadata[, cond := factor(cond, c('wild-type', 'knockout'))]
metadata = metadata[, .(sample, cond, time)]

idx = c(which(rownames(y) %in% c('13170', '12686', '26897')),
        round(seq(1, nrow(y), length.out = nGenes - 3)))
y = y[idx, metadata$sample]

qsave(metadata, file.path(dataDir, 'GSE34018_metadata.qs'))
qsave(y, file.path(dataDir, 'GSE34018_data.qs'))

# fit = getModelFit(y, metadata, nKnots = 3L, condColname = 'cond')
# fit = getPosteriorFit(fit)
# qsave(fit, file.path(dataDir, 'GSE34018_fit.qs'))
#
# rhyStats = getRhythmStats(fit)
# qsave(rhyStats, file.path(dataDir, 'GSE34018_rhy_stats.qs'))
#
# diffRhyStats = getDiffRhythmStats(fit, rhyStats, levels(metadata$cond))
# qsave(diffRhyStats, file.path(dataDir, 'GSE34018_diff_rhy_stats.qs'))
#
# fit = getPosteriorSamples(fit, nPosteriorSamples = 50)
# qsave(fit, file.path(dataDir, 'GSE34018_fit_ps.qs'))
