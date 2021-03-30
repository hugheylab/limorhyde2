#' @importFrom data.table data.table :=
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom zeallot %<-%
NULL


globalVariables(c(
  'condIdx', '.SD', 'peak_phase', 'trough_phase', 'rms_diff_rhy', 'nConds',
  'nK', 'nCon', 'nCov', 'lower', 'upper', 'nConds', 'co', 'peak_trough_amp',
  'rms_amp', 'mean_value', 'P.Value', 'adj.P.Val', '.', 'cond', 'feature',
  'nCovars', 'nKnots', 'peak_value', 'trough_value', 'period'))
  # 'ampl', 'cond', 'feature', 'nKnots', 'period', 'mNow', 'peak_value',
  # 'trough_value', 'cNum', 'condNow', 'peak_time', 'trough_time', '.'))


# addIntercept = function(b, intercept) {
#
#   stopifnot('intercept argument is either TRUE or FALSE'=is.logical(intercept))
#
#   if (intercept) {
#     b = cbind(1, b)
#     colnames(b)[1L] = 'intercept'}
#
#
#   return(b)}

addIntercept = function(b, intercept) {
  if (isTRUE(intercept)) {
    b = cbind(1, b)
    colnames(b)[1L] = 'intercept'}
  return(b)}


# getBasis = function(x, period = 24, nKnots = 4, intercept=FALSE){
#
#   stopifnot('Supply a numeric value > 0 for period' = length(period) == 1L,
#                                                         is.numeric(period),
#                                                         period > 0)
#
#   stopifnot('Supply a numeric value > 0 for nKnots' = length(nKnots) == 1L,
#             (is.numeric(nKnots) & nKnots >= 2) |is.null(nKnots))
#
#   if(is.null(nKnots)|identical(2, nKnots)){
#
#     b = cbind(cos(x / period * 2 * pi),
#               sin(x / period * 2 * pi))
#     colnames(b) = paste0('basis', 1:ncol(b))
#   } else {
#
#     knots = seq(0, period, length = nKnots + 2) # add 2 to include boundary values
#     b = pbs::pbs(x %% period, knots = knots[-c(1, length(knots))],
#                  Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
#     colnames(b) = paste0('basis', 1:nKnots)
#
#     b = b - (1/(nKnots+1)) }
#
#     b = addIntercept(b, intercept)
#
#
#     return(b)}

#' @export
getBasis = function(time, period = 24, nKnots = 4, intercept = TRUE) {
  stopifnot(length(period) == 1L,
            is.numeric(period),
            period > 0,
            length(nKnots) == 1L,
            is.numeric(nKnots),
            is.null(nKnots) || nKnots >= 2,
            length(intercept) == 1L,
            is.logical(intercept))

  if (is.null(nKnots)) nKnots = 2

  if (nKnots == 2) {
    tt = time / period * 2 * pi
    b = cbind(cos(tt), sin(tt))
  } else {
    knots = seq(0, period, length = nKnots + 2) # including boundary knots
    b = pbs::pbs(time %% period, knots = knots[-c(1, length(knots))],
                 Boundary.knots = knots[c(1, length(knots))])[, , drop = FALSE]
    # scale basis so intercept doesn't change when other coefs shrink
    b = b - 1 / (nKnots + 1)}

  colnames(b) = paste0('basis', 1:nKnots)
  b = addIntercept(b, intercept)
  return(b)}


# getSm = function(md, timeColname, conditionsColname){
#
#   md = as.data.table(md)
#   stopifnot( 'time column is not in metadata'= length(timeColname) == 1L,
#              'time column is not in metadata'= timeColname %in% colnames(md),
#              is.character(timeColname))
#
#   if(is.null(conditionsColname)){
#
#     colsKeep = c(timeColname)
#     sm = md[, colsKeep, with=FALSE]
#     setnames(sm, colsKeep, 'time')
#     sm[, cond := 'one']} else{
#
#       stopifnot('Specify a string indicating the condition/treatment column name in metadata'=
#                   length(conditionsColname) == 1L,
#                 'A column indicating the condition for each sample is not in metadata'=
#                   conditionsColname %in% colnames(md),
#                 is.character(conditionsColname))
#       colsKeep = c(timeColname, conditionsColname)
#
#       sm = md[, colsKeep, with=FALSE] # or use :=
#       setnames(sm, colsKeep, c('time','cond'))
#   }
#
#   stopifnot('time column must have numbers only' = is.numeric(sm$time))
#
#   return(sm)
# }

getMetadata = function(metadata, timeColname, condColname, covarColnames) {
  stopifnot(length(timeColname) == 1L,
            timeColname %in% colnames(metadata),
            is.numeric(metadata[[timeColname]]),
            is.null(condColname) || length(condColname) == 1L)

  m = data.table(time = metadata[[timeColname]])

  if (!is.null(condColname)) {
    stopifnot(condColname != timeColname,
              condColname %in% colnames(metadata))
    m[, cond := factor(metadata[[condColname]])]}

  if (!is.null(covarColnames)) {
    covarColnames = unique(covarColnames)
    stopifnot(!any(covarColnames %in% c(timeColname, condColname)),
              all(covarColnames %in% colnames(metadata)))
    for (covarName in covarColnames) {
      data.table::set(
        m, NULL, paste0('covar_', covarName), metadata[[covarName]])}}

  return(m[])}


getDesign = function(metadata, period, nKnots) {
  b = getBasis(metadata$time, period, nKnots, FALSE)
  m = cbind(metadata, b)

  r = paste0('basis', 1:nKnots, collapse = ' + ')
  if ('cond' %in% colnames(m)) r = sprintf('cond * (%s)', r)

  covarIdx = startsWith(colnames(m), 'covar')
  if (any(covarIdx)) r = paste(c(r, colnames(m)[covarIdx]), collapse = ' + ')
  design = stats::model.matrix(stats::formula(paste('~', r)), data = m)

  covarIdx = startsWith(colnames(design), 'covar')
  if (any(covarIdx)) design = design[, c(which(!covarIdx), which(covarIdx))]

  nConds = length(unique(m$cond)) # works even if cond is not a column
  if (nConds > 2) {
    nCols = ncol(design) - sum(covarIdx)
    idx = foreach(condIdx = 2:nConds, .combine = c) %do% {
      seq(c(nConds + nKnots + condIdx - 1), nCols, nConds - 1)}
    idx = c(1:(nConds + nKnots), idx)
    if (any(covarIdx)) idx = c(idx, (nCols + 1):ncol(design))
    design = design[, idx]}
  return(design)}


# getCK = function(mat){
#
#   cols = colnames(mat)
#   nCond = which(cols == 'basis1') - 1
#   nKnots = (length(cols) - nCond)/nCond
#
#   ck = c(nCond, nKnots)
#
#   return(ck)
#
# }

getNumKnotCondCovar = function(cols) {
  nCovars = sum(startsWith(cols, 'covar'))
  nConds = which(cols == 'basis1') - 1L
  nKnots = as.integer((length(cols) - nConds - nCovars) / nConds)
  return(c(nKnots, nConds, nCovars))}
