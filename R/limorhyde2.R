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


getCosinorBasisAsh = function(fit, mixcompdist = 'halfnormal', ci = FALSE, ...){
  #get modified cost and sint coefs for a single condition
  co = fit$coefficients

  i = 1L
  c(ix, iy) %<-% fit$coef_lookup[cond_num == i  & model_comp != 'intercept']$coef_idx

  # statsDT[, sumCost := rowSums(.SD), .SDcols = cosIdx]
  # statsDT[, sumSint := rowSums(.SD), .SDcols = sinIdx]

  statsDT = data.table(fit$coefficients, keep.rownames = "feature")
  statsDT = melt(statsDT, measure = c("cost", "sint"), value.name = "estimate")

  formList = list(~ x1, ~x2)
  idx = c(ix, iy)

  seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                     fit$cov.coefficients[idx, idx])

  statsDT[, se := as.vector(seMat)]

  ashObj = statsDT[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
  d = cbind(statsDT[, .(feature, variable)], setDT(ashObj$result))

  return(d)}




getCosinorBasisAsh2 = function(fit, mixcompdist = 'halfnormal', ci = FALSE, ...){
  #get modified cost and sint coefs for a multiple conditions
  co = fit$coefficients
  statsDT = data.table(fit$coefficients, keep.rownames = "feature")

  #get modified coefficients for k = 0
  i = 1L
  c(ix, iy) %<-% fit$coef_lookup[cond_num == i  & model_comp != 'intercept']$coef_idx

  formList = list(~ x1, ~x2)
  idx = c(ix, iy)

  seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                     fit$cov.coefficients[idx, idx])
  d0 = data.table(cond = rep( fit$coef_lookup[cond_num == i]$cond_lev[1L],
                              length.out = nrow(co)))
  d0 = melt(statsDT[,.(feature, cost, sint)],
            measure = c("cost", "sint"),
            value.name = "estimate")
  d0[,cond := fit$coef_lookup[cond_num == i]$cond_lev[1L]]
  d0[, se := as.vector(seMat)]
  ashObj = d0[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
  d0 = cbind(d0, setDT(ashObj$result))



  # get modified coefficients for k = 1
  c(ix, iy) %<-% fit$coef_lookup[model_comp %like% 'cost']$coef_idx
  c(jx, jy) %<-% fit$coef_lookup[model_comp %like% 'sint']$coef_idx

  formList = list(~ x1 + x2, ~x3 + x4)
  idx = c(ix, iy, jx, jy)

  seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                     fit$cov.coefficients[idx, idx])

  statsDT[, sum_cost := co[, ix] + co[, iy]]
  statsDT[, sum_sint := co[jx] + co[, jy]]
  d1 = melt(statsDT[, .(feature, sum_cost, sum_sint)],
            measure = c("sum_cost", "sum_sint"),
            value.name = "estimate")
  d1[, cond := fit$coef_lookup[!(cond_num == i)]$cond_lev[1L]]
  d1[, se := as.vector(seMat)]

  ashObj = d1[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
  d1 = cbind(d1, setDT(ashObj$result))

  d = rbind(d0,d1)

  return(d)}

#' @export
getRhythmStats = function(fit) {
  co = fit$coefficients
  # levs = fit$cond_levels

  i = 1L
  d1 = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[1L],# levs[i],
                  feature = rep(rownames(co), 3),
                  variable = rep(c('mesor', 'amp', 'phase'), each = nrow(co)))

  c(ix, iy) %<-% fit$coef_lookup[cond_num == i & model_comp != 'intercept']$coef_idx # fit$time_idx
  d1[variable == 'mesor', estimate := co[, i]]
  d1[variable == 'amp', estimate := sqrt(co[, ix]^2 + co[, iy]^2)]
  d1[variable == 'phase', estimate := atan2(co[, iy], co[, ix])]

  formList = list(~ x1, ~ sqrt(x2^2 + x3^2), ~ atan2(x3, x2))
  idx = c(i, ix, iy)

  seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                     fit$cov.coefficients[idx, idx])
  seMat[, 3] = seMat[, 3] * fit$period / 2 / pi
  d1[, se := as.vector(seMat)]

  d2 = foreach(j = 2:length(levels(fit$coef_lookup$cond_lev)),
               .combine = rbind) %dopar% {

    dj = data.table(cond = fit$coef_lookup[cond_num == j]$cond_lev[1L],
                    feature = rep(rownames(co), 3),
                    variable = rep(c('mesor', 'amp', 'phase'), each = nrow(co)))

    # c(jx, jy) %<-% fit$cond_time_idx[(2 * (j - 1) - 1):(2 * (j - 1))]
    c(jx, jy) %<-% fit$coef_lookup[cond_num == j & model_comp != 'intercept']$coef_idx

    dj[variable == 'mesor',
       estimate := co[, i] + co[, j]]
    dj[variable == 'amp',
       estimate := sqrt((co[, ix] + co[, jx])^2 + (co[, iy] + co[, jy])^2)]
    dj[variable == 'phase',
       estimate := atan2(co[, iy] + co[, jy], co[, ix] + co[, jx])]

    formList = list(~ x1 + x2,
                    ~ sqrt((x3 + x5)^2 + (x4 + x6)^2),
                    ~ atan2(x4 + x6, x3 + x5))
    idx = c(i, j, ix, iy, jx, jy)

    seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                       fit$cov.coefficients[idx, idx])
    seMat[, 3] = seMat[, 3] * fit$period / 2 / pi
    dj[, se := as.vector(seMat)]}

  d3 = rbind(d1, d2)
  d3[variable == 'phase' & estimate < 0, estimate := estimate + 2 * pi]
  d3[variable == 'phase' & estimate >= 2 * pi, estimate := estimate - 2 * pi]
  d3[variable == 'phase', estimate := estimate * fit$period / 2 / pi]
  return(d3[])}




getRhythmStats2 = function(fit, basisAsh){

  i = 1
  coPost = dcast(basisAsh[, .(feature, variable, PosteriorMean, PosteriorSD)],
                 feature ~ variable, value.var = c("PosteriorMean"))

  coPost[, intercept := fit$coefficients[,1L]]
  coPost[, condwt := fit$coefficients[,2L]]

  coMat = as.matrix(coPost, rownames = "feature")

  # calculate amplitude and phase for k = 0
  d0 = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[i],
                  feature = coPost$feature,
                  variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coPost)))

  d0[variable == 'mesor', estimate := coMat[, 5]]
  d0[variable == 'amp', estimate := sqrt(coMat[, 1]^2 + coMat[, 2]^2)]
  d0[variable == 'phase', estimate := atan2(coMat[, 2], coMat[, 1])]

  formList = list(~ x1, ~ sqrt(x2^2 + x3^2), ~ atan2(x3, x2))

  #calculate amplitude and phase for k = 1

  d1 = data.table(cond = fit$coef_lookup[!(cond_num == i)]$cond_lev[i],
                  feature = coPost$feature,
                  variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coPost)))

  d1[variable == 'mesor', estimate := coMat[, 6]]
  d1[variable == 'amp', estimate := sqrt(coMat[, 3]^2 + coMat[, 4]^2)]
  d1[variable == 'phase', estimate := atan2(coMat[, 4], coMat[, 3])]



  #calculate amplitude and phase difference
  d = rbind(d0, d1)

  return(d[])}




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




getOneCondStats = function(fit, basisAsh, confLevel = 0.95){
  #single condition stats
  co = fit$coefficients
  i = 1
  c(ix, iy) %<-% fit$coef_lookup[cond_num == i & model_comp != 'intercept']$coef_idx # fit$time_idx
  idx = c(i, ix, iy)

  coPost = dcast(basisAsh[, .(feature, variable, PosteriorMean, PosteriorSD)],
                 feature ~ variable, value.var = c("PosteriorMean"))
  coPost[, intercept := fit$coefficients[,i]]
  coPost = coPost[, c(1,4,2,3)] #rearrange cols to fit indx

  coMat = as.matrix(coPost, rownames = "feature")

  d1 = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[i],
                  feature = coPost$feature,
                  variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coPost)))

  d1[variable == 'mesor', estimate := coMat[, i]]
  d1[variable == 'amp', estimate := sqrt(coMat[, ix]^2 + coMat[, iy]^2)]
  d1[variable == 'phase', estimate := atan2(coMat[, iy], coMat[, ix])]

  #use delta method on se vals from ashr
  formList = list(~ x1, ~ sqrt(x2^2 + x3^2), ~ atan2(x3, x2))

  seMat = getStdErrs(formList, coMat, fit$s2.post,
                     fit$cov.coefficients[idx, idx])
  seMat[, 3] = seMat[, 3] * fit$period / 2 / pi
  d1[, se := as.vector(seMat)]


  d1[variable == 'phase' & estimate < 0, estimate := estimate + 2 * pi]
  d1[variable == 'phase' & estimate >= 2 * pi, estimate := estimate - 2 * pi]
  d1[variable == 'phase', estimate := estimate * fit$period / 2 / pi]


  quant = stats::qt(0.5 + confLevel / 2, fit$df.total)
  d1[, `:=`(ci_lower = estimate - se * quant,
            ci_upper = estimate + se * quant,
            tstat = estimate / se,
            pval = 2 * stats::pt(-abs(estimate / se), df = rep(fit$df.total, 2)))]


  return(d1)}


#' @export
getDiffRhythmStats = function(fit, rhythmStats,
                              conds = levels(fit$coef_lookup$cond_lev)[1:2],
                              confLevel = 0.95) {
  stopifnot(length(conds) == 2, all(conds %in% levels(fit$coef_lookup$cond_lev)),
            length(confLevel) == 1, is.numeric(confLevel),
            confLevel > 0, confLevel < 1)

  d0 = rhythmStats[cond %in% conds]
  d0[, cond := factor(cond, conds)]
  setorderv(d0, 'cond')

  d1 = d0[, .(estimate = diff(estimate)), by = .(feature, variable)]
  d1[, variable := paste0('diff_', variable)]

  d1[variable == 'diff_phase', estimate := estimate %% fit$period]
  d1[variable == 'diff_phase' & estimate > fit$period / 2,
     estimate := estimate - fit$period]
  d1[variable == 'diff_phase' & estimate < (-fit$period / 2),
     estimate := estimate + fit$period]

  levs = levels(fit$coef_lookup$cond_lev)

  if (conds[1L] == levs[1L]) { #fit$cond_levels[1]) {
    # cos1 -> x2
    # sin1 -> x3
    # cos2 -> x2 + x4
    # sin2 -> x3 + x5
    formList = list(~ x1,
                    ~ sqrt((x2 + x4)^2 + (x3 + x5)^2) - sqrt(x2^2 + x3^2),
                    ~ atan2(x2 * (x3 + x5) - x3 * (x2 + x4),
                            x2 * (x2 + x4) + x3 * (x3 + x5)))

    # idxTmp = match(conds[2], fit$cond_levels) - 1
    # idx = c(idxTmp + 1, fit$time_idx,
    #         fit$cond_time_idx[(2 * idxTmp - 1):(2 * idxTmp)])
    idx = c(fit$coef_lookup[cond_lev == conds[2L] & model_comp == 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[1L] & model_comp != 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[2L] & model_comp != 'intercept']$coef_idx)

  } else if (conds[2L] == levs[1L]) { #fit$cond_levels[1]) {
    # cos1 -> x2 + x4
    # sin1 -> x3 + x5
    # cos2 -> x2
    # sin2 -> x3
    formList = list(~ -x1,
                    ~ sqrt(x2^2 + x3^2) - sqrt((x2 + x4)^2 + (x3 + x5)^2),
                    ~ atan2((x2 + x4) * x3 - (x3 + x5) * x2,
                            x2 * (x2 + x4) + x3 * (x3 + x5)))

    # idxTmp = match(conds[2], fit$cond_levels) - 1
    # idx = c(idxTmp + 1, fit$time_idx,
    #         fit$cond_time_idx[(2 * idxTmp - 1):(2 * idxTmp)])
    idx = c(fit$coef_lookup[cond_lev == conds[2L] & model_comp == 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[1L] & model_comp != 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[2L] & model_comp != 'intercept']$coef_idx)

  } else {
    # cos1 -> x3 + x5
    # sin1 -> x4 + x6
    # cos2 -> x3 + x7
    # sin2 -> x4 + x8
    formList = list(~ x2 - x1,
                    ~ sqrt((x3 + x5)^2 + (x4 + x6)^2) -
                      sqrt((x3 + x7)^2 + (x4 + x8)^2),
                    ~ atan2(x3 + x5 * (x4 + x8) - (x4 + x6) * (x3 + x7),
                            (x3 + x5) * (x3 + x7) + (x4 + x6) * (x4 + x8)))

    # idxTmp = match(conds, fit$cond_levels) - 1
    # idx = c(idxTmp + 1, fit$time_idx,
    #         fit$cond_time_idx[(2 * idxTmp[1] - 1):(2 * idxTmp[1])],
    #         fit$cond_time_idx[(2 * idxTmp[2] - 1):(2 * idxTmp[2])])
    idx = c(fit$coef_lookup[cond_lev == conds[1L] & model_comp == 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[2L] & model_comp == 'intercept']$coef_idx,
            fit$coef_lookup[cond_num == 1L & model_comp != 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[1L] & model_comp != 'intercept']$coef_idx,
            fit$coef_lookup[cond_lev == conds[2L] & model_comp != 'intercept']$coef_idx)}

  seMat = getStdErrs(formList, fit$coefficients[, idx],
                     fit$s2.post, fit$cov.coefficients[idx, idx])
  seMat[, 3L] = seMat[, 3L] * fit$period / 2 / pi
  d2 = data.table(feature = rep(rownames(fit), 3),
                  variable = rep(c('diff_mesor', 'diff_amp', 'diff_phase'), each = nrow(fit)),
                  se = as.vector(seMat))

  quant = stats::qt(0.5 + confLevel / 2, fit$df.total)
  d = merge(d1, d2, by = c('feature', 'variable'), sort = FALSE)
  d[, `:=`(ci_lower = estimate - se * quant,
           ci_upper = estimate + se * quant,
           tstat = estimate / se,
           pval = 2 * stats::pt(-abs(estimate / se), df = rep(fit$df.total, 2)))]
  return(d[])}



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
getDiffRhythmTests = function(fit, diffRhythmStats = NULL, qvalRhyCutoff = 0.2) {
  topTables = getTopTables(fit, qvalRhyCutoff)
  d1 = merge(topTables$cond_time[, .(feature, pval_diff_rhy = P.Value, qval_diff_rhy = adj.P.Val)],
             topTables$time[, .(feature, pval_rhy = P.Value, qval_rhy = adj.P.Val)],
             by = 'feature', sort = FALSE)

  if (is.null(diffRhythmStats)) {
    setorderv(d1, c('pval_diff_rhy', 'pval_rhy'))
    return(d1[])}

  d2 = data.table::dcast(diffRhythmStats, feature ~ variable, value.var = 'pval')
  setnames(d2, function(x) paste0('pval_', x))
  setnames(d2, 'pval_feature', 'feature')

  d = merge(d1, d2, by = 'feature')
  d[, qval_diff_mesor := stats::p.adjust(pval_diff_mesor, 'BH')]
  d[!is.na(qval_diff_rhy), qval_diff_amp := stats::p.adjust(pval_diff_amp, 'BH')]
  d[!is.na(qval_diff_rhy), qval_diff_phase := stats::p.adjust(pval_diff_phase, 'BH')]

  setcolorder(d, c('feature', 'pval_diff_rhy', 'qval_diff_rhy',
                   'pval_diff_amp', 'qval_diff_amp', 'pval_diff_phase', 'qval_diff_phase',
                   'pval_rhy', 'qval_rhy', 'pval_diff_mesor', 'qval_diff_mesor'))
  setorderv(d, c('pval_diff_rhy', 'pval_rhy', 'pval_diff_mesor'))
  return(d[])}


#' @export
getDiffRhythmAsh = function(diffRhythmStats, diffRhythmTests,
                            mixcompdist = 'halfnormal', ci = FALSE, ...) {

  # mesor
  dNow = diffRhythmStats[variable == 'diff_mesor']
  ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]
  dm = cbind(dNow[, .(feature, variable)], setDT(ashObj$result))
  if (isTRUE(ci)) {
    ciMat = ashr::ashci(ashObj)
    dm[, `:=`(credint_lower = ciMat[, 1L], credint_upper = ciMat[, 2L])]}

  dNow = merge(diffRhythmStats[variable != 'diff_mesor',
                               .(feature, variable, estimate, se)],
               diffRhythmTests[!is.na(qval_diff_rhy), .(feature)],
               by = 'feature')

  # amplitude
  ashObj = dNow[variable == 'diff_amp',
                ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]

  da = cbind(dNow[variable == 'diff_amp', .(feature, variable)],
             setDT(ashObj$result))
  if (isTRUE(ci)) {
    ciMat = ashr::ashci(ashObj)
    da[, `:=`(credint_lower = ciMat[, 1L], credint_upper = ciMat[, 2L])]}

  # phase
  ashObj = dNow[variable == 'diff_phase',
                ashr::ash(estimate, se, mixcompdist = mixcompdist, ...)]

  dp = cbind(dNow[variable == 'diff_phase', .(feature, variable)],
             setDT(ashObj$result))
  if (isTRUE(ci)) {
    ciMat = ashr::ashci(ashObj)
    dp[, `:=`(credint_lower = ciMat[, 1L], credint_upper = ciMat[, 2L])]}

  d = rbind(dm, da, dp)
  return(d)}

