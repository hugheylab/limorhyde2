library('data.table')
library('qs')

rawDir = 'data-raw'
nGenes = 50L

########################################

metadata = fread(file.path(rawDir, 'GSE54650_sample_metadata.csv'))
y = qread(file.path(rawDir, 'GSE54650_expression_data.qs'))

metadata = metadata[tissue == 'liver', .(sample, tissue, time)]

idx = c(which(rownames(y) %in% c('13088', '13170', '13869')),
        round(seq(1, nrow(y), length.out = nGenes - 3)))
y = y[idx, metadata$sample]

GSE54650 = list(y = y, metadata = metadata)

usethis::use_data(GSE54650, overwrite = TRUE)

########################################

metadata = fread(file.path(rawDir, 'GSE34018_sample_metadata.csv'))
y = qread(file.path(rawDir, 'GSE34018_expression_data.qs'))

metadata[, cond := factor(cond, c('wild-type', 'knockout'))]
metadata = metadata[, .(sample, cond, time)]

idx = c(which(rownames(y) %in% c('13170', '12686', '26897')),
        round(seq(1, nrow(y), length.out = nGenes - 3)))
y = y[idx, metadata$sample]

GSE34018 = list(y = y, metadata = metadata)

usethis::use_data(GSE34018, overwrite = TRUE)
