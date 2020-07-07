#' @importFrom data.table data.table := setnames setorderv setcolorder setDT
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom zeallot %<-%
NULL


globalVariables(c('estimate', 'feature', 'variable', 'se', 'meanNow', 'j', 'ix',
                  'iy', 'jx', 'jy', 'pval_diff_amp', 'pval_diff_mesor',
                  'pval_diff_phase', 'P.Value', 'adj.P.Val', 'cond', 'varNow',
                  'qval_diff_phase', 'qval_diff_rhy', 'pval_diff_rhy',
                  '.', 'pval_rhy', 'qval_diff_amp', 'qval_diff_mesor',
                  'cond_lev', 'cond_num', 'model_comp', 'coef_idx'))


addIntercept = function(b, intercept) {
  if (intercept) {
    b = cbind(1, b)
    colnames(b)[1L] = 'intercept'}
  return(b)}


getCosinorBasis = function(x, period, intercept) {
  b = cbind(cos(x / period * 2 * pi),
            sin(x / period * 2 * pi))
  colnames(b) = c('cos', 'sin')
  b = addIntercept(b, intercept)
  return(b)}


#' @export

getModelFit = function(x, metadata, condColname = NULL, timeColname = 'time',
                       period = 24, formSupp = NULL, ...) {
  # columns of x must correspond to rows of metadata
  stopifnot(length(period) == 1L, is.numeric(period), period > 0)

  metadata = as.data.table(metadata)

  if (!is.null(condColname)) {
    setnames(metadata, c(condColname, timeColname), c('cond', 'time'))
    if (!is.factor(metadata$cond)) {
      metadata[, cond := factor(cond)]}
  } else {
    # metadata = metadata[, .()]  #make sure time and sample_id columns are present
    setnames(metadata, timeColname, 'time')}


  b = getCosinorBasis(metadata$time, period, FALSE)
  colnames(b) = paste0(colnames(b), 't')
  metadata = cbind(metadata, b)

  if (is.null(formSupp)) {

    if (is.null(condColname)){
      formFull = ~ cost + sint
    } else{
      formFull = ~ cond * (cost + sint)}

  }else {
    # must be a formula that has no interactions and no mention of cond or time
    stopifnot(inherits(formSupp, 'formula'),
              !any(grepl(':', attr(stats::terms(formSupp), 'term.labels'))),
              !any(attr(stats::terms(formSupp), 'term.labels') %in%
                     c(condColname, timeColname)))
    formFull = stats::update.formula(formSupp, ~ . + cond * (cost + sint))}

  design = stats::model.matrix(formFull, data = metadata)
  if (!is.null(formSupp)) {
    # won't work if formSupp has interactions
    nSupp = length(attr(stats::terms(formSupp), 'term.labels'))
    idxNew = c(1L, (2L + nSupp):ncol(design), 2:(1L + nSupp))
    design = design[, idxNew]}

  fit = limma::lmFit(x, design)
  fit = limma::eBayes(fit, trend = TRUE, ...)

  if(is.null(condColname)){ levs = 1L
  }else{ levs = levels(metadata$cond)}
  fit$coef_lookup = data.table(cond_lev = rep(levs, each = 3),
                               model_comp = rep(c('intercept', 'cost', 'sint'),
                                                length(levs)))
  fit$coef_lookup[, cond_lev := factor(cond_lev, levs)]
  fit$coef_lookup[, cond_num := as.numeric(cond_lev)]
  fit$coef_lookup[cond_num == 1L, coef_idx := c(1, length(levs) + 1:2)]
  if (!is.null(condColname)){
    for (i in 2:length(levs)) {
      fit$coef_lookup[cond_num == i,
                      coef_idx := c(i, length(levs) + i + c(1, length(levs)))]}}


  fit$period = period
  return(fit)}



getTopTables = function(fit, qvalRhyCutoff = 0.2, sort = TRUE) {
  coefIdx = fit$coef_lookup[model_comp == 'intercept' & cond_num != 1L]$coef_idx
  dCond = limma::topTable(fit, coef = coefIdx, number = Inf, sort.by = 'none',
                          confint = TRUE)
  # confint only returned if <= 2 levels
  dCond = data.table(dCond, keep.rownames = TRUE)
  setnames(dCond, 'rn', 'feature')

  coefIdx = fit$coef_lookup[model_comp != 'intercept']$coef_idx
  dTime = limma::topTable(fit, coef = coefIdx, number = Inf, sort.by = 'none')
  dTime = data.table(dTime, keep.rownames = TRUE)
  setnames(dTime, 'rn', 'feature')

  coefIdx = fit$coef_lookup[model_comp != 'intercept' & cond_num != 1L]$coef_idx
  dCondTime = limma::topTable(fit, coef = coefIdx, number = Inf, sort.by = 'none')
  dCondTime = data.table(dCondTime, keep.rownames = TRUE)
  setnames(dCondTime, 'rn', 'feature')

  dCondTime[, adj.P.Val := NA_real_]
  dCondTime[dTime[adj.P.Val <= qvalRhyCutoff, which = TRUE],
            adj.P.Val := stats::p.adjust(P.Value, 'BH')]

  if (isTRUE(sort)) {
    setorderv(dCond, 'P.Value')
    setorderv(dTime, 'P.Value')
    setorderv(dCondTime, 'P.Value')}

  return(list(cond = dCond, time = dTime, cond_time = dCondTime))}


getStdErrs = function(g, meanEst, varEst, covBase) {
  syms = paste0('x', 1:ncol(meanEst))
  gDeriv = lapply(g, function(form) Deriv::Deriv(form, syms))
  se = foreach(meanNow = iterators::iter(meanEst, by = 'row'), varNow = varEst,
               .combine = rbind) %dopar% {
    for (i in 1:length(meanNow)) {
      assign(syms[i], meanNow[i])}
    gdashmu = sapply(gDeriv, function(e) eval(e))
    covNow = t(gdashmu) %*% (varNow * covBase) %*% gdashmu
    seNow = sqrt(diag(covNow))}
  return(se)}



#get ashr modified estimates for sint and cost

getCosinorBasisAsh2 = function(fit, mixcompdist = 'halfnormal', ci = FALSE, ...){

  co = fit$coefficients
  statsDT = data.table(fit$coefficients, keep.rownames = "feature")

  #get modified coefficients for k = 0
  i = 1L
  c(ix, iy) %<-% fit$coef_lookup[cond_num == i  & model_comp != 'intercept']$coef_idx

  formList = list(~ x1, ~x2)
  idx = c(ix, iy)

  seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                     fit$cov.coefficients[idx, idx])
  d0 = melt(statsDT[,.(feature, cost, sint)],
            measure = c("cost", "sint"),
            value.name = "estimate")
  d0[,cond := fit$coef_lookup[cond_num == i]$cond_lev[i]]
  d0[, se := as.vector(seMat)]
  ashObj = d0[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
  d0 = cbind(d0, setDT(ashObj$result))


  # get modified coefficients for k = 1
  if (length(levels(fit$coef_lookup$cond_lev)) > 1){
    d1 = foreach(j = 2:length(levels(fit$coef_lookup$cond_lev)),
                 .combine = rbind) %dopar% {

                   c(ix, iy) %<-% fit$coef_lookup[cond_num %in% c(i,j) & model_comp %like% 'cost']$coef_idx
                   c(jx, jy) %<-% fit$coef_lookup[cond_num %in% c(i,j) & model_comp %like% 'sint']$coef_idx

                   formList = list(~ x1 + x2, ~x3 + x4)
                   idx = c(ix, iy, jx, jy)

                   seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                                      fit$cov.coefficients[idx, idx])

                   statsDT[, cost_k1 := co[, ix] + co[, iy]]
                   statsDT[, sint_k1 := co[jx] + co[, jy]]
                   dNow = melt(statsDT[, .(feature, cost = cost_k1, sint = sint_k1)],
                               measure = c("cost", "sint"),
                               value.name = "estimate")
                   dNow[, cond := fit$coef_lookup[cond_num == j]$cond_lev[1L]]
                   dNow[, se := as.vector(seMat)]

                   ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
                   dRes = cbind(dNow, setDT(ashObj$result))}

    d = rbind(d0,d1)} else{

      d = d0}

  return(d)}


getRhythmStats2 = function(fit, basisAsh){

  i = 1L
  coRaw = fit$coefficients

  # calculate amplitude and phase for k = 0
  c(ix, iy, iz) %<-% fit$coef_lookup[cond_num == i]$coef_idx
  dRaw = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[i],
                    feature = rep(rownames(coRaw), 3),
                    variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coRaw)))
  #how to calculate mesor from sin and cos coefs?
  dRaw[variable == 'mesor', estimate_raw := coRaw[,ix]]
  dRaw[variable == 'amp', estimate_raw := sqrt(coRaw[, iy]^2 + coRaw[, iz]^2)]
  dRaw[variable == 'phase', estimate_raw := atan2(coRaw[, iz], coRaw[, iy])]

  # calculate moderated values for rhythmic parameters
  coPost = dcast(basisAsh[, .(feature, cond, variable, PosteriorMean)],
                 feature ~ variable + cond, value.var = c("PosteriorMean"))

  dPost = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[i],
                     feature = rep(rownames(coRaw), 3),
                     variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coRaw)))

  dPost[variable == 'mesor', estimate_mod := coRaw[,i]]
  dPost[variable == 'amp', estimate_mod := sqrt(coPost[, i+1]^2 + coPost[, i+2]^2)]
  dPost[variable == 'phase', estimate_mod := atan2(coPost[, i+2], coPost[, i+1])]

  d0 = merge(dRaw, dPost, by = c("feature", "variable", "cond"))

  #calculate amplitude and phase for k = 1
  if (length(levels(fit$coef_lookup$cond_lev)) >= 2) {
    d1 = foreach(j = 2:length(levels(fit$coef_lookup$cond_lev)),
                 .combine = rbind) %do% {

                   c(jx, jy, jz) %<-% fit$coef_lookup[cond_num == j]$coef_idx
                   dRaw = data.table(cond = fit$coef_lookup[cond_num == j]$cond_lev[i],
                                     feature = rep(rownames(coRaw), 3),
                                     variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coPost)))

                   dRaw[variable == 'mesor',
                        estimate_raw := coRaw[, ix] + coRaw[, jx]]
                   dRaw[variable == 'amp',
                        estimate_raw := sqrt((coRaw[, iy] + coRaw[, jy])^2 + (coRaw[, iz] + coRaw[jz])^2)]
                   dRaw[variable == 'phase',
                        estimate_raw := atan2(coRaw[, iz] + coRaw[, jz], coRaw[, iy] + coRaw[jy])]


                   x = colnames(coPost) %like% "feature" |
                     colnames(coPost) %like% fit$coef_lookup[cond_num == j]$cond_lev[i]
                   coPost = coPost[,(x), with = FALSE]
                   coPost = as.matrix(coPost, rownames = 'feature')

                   dPost = data.table(cond = fit$coef_lookup[cond_num == j]$cond_lev[i],
                                      feature = rep(rownames(coRaw), 3),
                                      variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coRaw)))

                   dPost[variable == 'mesor',
                         estimate_mod := coRaw[, ix] + coRaw[, jx]]
                   dPost[variable == 'amp',
                         estimate_mod := sqrt(coPost[, 1]^2 + coPost[, 2]^2)]
                   dPost[variable == 'phase', estimate_mod := atan2(coPost[, 2], coPost[, 1])]

                   tl = merge(dRaw, dPost, by = c("feature", "variable", "cond"))}

    #calculate amplitude and phase difference
    d = rbind(d0, d1)}

  d[variable == 'phase' & estimate_raw < 0, estimate_raw := estimate_raw + 2 * pi]
  d[variable == 'phase' & estimate_raw >= 2 * pi, estimate_raw := estimate_raw - 2 * pi]
  d[variable == 'phase', estimate_raw := estimate_raw * fit$period / 2 / pi]


  d[variable == 'phase' & estimate_mod < 0, estimate_mod := estimate_mod + 2 * pi]
  d[variable == 'phase' & estimate_mod >= 2 * pi, estimate_mod := estimate_mod - 2 * pi]
  d[variable == 'phase', estimate_mod := estimate_mod * fit$period / 2 / pi]


  return(d[])}



#' @export

getDiffRhythmStats2 =  function(fit, rhythmStats,
                                conds = levels(fit$coef_lookup$cond_lev)[1:2],
                                confLevel = 0.95){

  stopifnot(length(conds) == 2, all(conds %in% levels(fit$coef_lookup$cond_lev)),
            length(confLevel) == 1, is.numeric(confLevel),
            confLevel > 0, confLevel < 1)

  d1 = rhythmStats[, .(estimate = diff(estimate)), by = .(feature, variable)]


  d1[variable == 'diff_phase', estimate := estimate %% fit$period]
  d1[variable == 'diff_phase' & estimate > fit$period / 2,
     estimate := estimate - fit$period]
  d1[variable == 'diff_phase' & estimate < (-fit$period / 2),
     estimate := estimate + fit$period]

  return(d1)}



#' @export

getRhythmTests = function(fit, oneCondStats, qvalRhyCutoff = 0.2){
  # one condition
  coefIdx = fit$coef_lookup[model_comp != 'intercept']$coef_idx
  dTime = limma::topTable(fit, coef = coefIdx, number = Inf, sort.by = 'none')
  dTime = data.table(dTime, keep.rownames = TRUE)
  setnames(dTime, 'rn', 'feature')

  d1 = dTime[, .(feature, pval_rhy = P.Value, qval_rhy = adj.P.Val)]

  # oneCondStats[, qval := stats::p.adjust(pval, 'BH')]
  d2 = oneCondStats[, qval := stats::p.adjust(pval, 'BH')]
  d3 = data.table::dcast(d2, feature ~ variable, value.var = c('pval', 'qval'))

  d = merge(d1, d3, by = 'feature')
  return(d)}



