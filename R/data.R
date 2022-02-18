#' Gene expression data for GSE54650
#'
#' Data are based on total RNA, measured by microarray, obtained from livers of
#' wild-type mice at various times after transfer to constant darkness. To save
#' space and time, the data include only a subset of genes, and so are mainly
#' useful for examples of how to use `limorhyde2`.
#'
#' @format A list with two elements:
#'
#' * `y`: Matrix of normalized, log-transformed expression values. Rows
#'   correspond to genes (rownames are Entrez Gene IDs) and columns to samples.
#' * `metadata`: data.table with one row per sample. `time` is in hours.
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54650>
#'
#' @seealso [GSE34018], [getModelFit()]
'GSE54650'


#' Gene expression data for GSE34018
#'
#' Data are based on total RNA, measured by microarray, obtained from livers of
#' wild-type and liver-specific Reverb alpha/beta double knockout mice at
#' various times in a 12h:12h light:dark cycle. To save space and time, the data
#' include only a subset of genes, and so are mainly useful for examples of how
#' to use `limorhyde2`.
#'
#' @format A list with two elements:
#'
#' * `y`: Matrix of normalized, log-transformed expression values. Rows
#'   correspond to genes (rownames are Entrez Gene IDs) and columns to samples.
#' * `metadata`: data.table with one row per sample. `time` is in hours.
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34018>
#'
#' @seealso [GSE54650], [getModelFit()]
'GSE34018'
