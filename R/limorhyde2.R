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

getCosinorBasisAsh2 = function(fit, combine_dist = TRUE, mixcompdist = 'halfnormal', ci = FALSE, ...){

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
  # d0[, se :=  rep(fit$s2.post,2)]

  if(isTRUE(combine_dist)){
    ashObj = d0[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
    d1 = cbind(d0, setDT(ashObj$result)) } else {

      dNow = d0[variable == 'cost']
      ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
      d_cost = cbind(dNow, setDT(ashObj$result))


      dNow = d0[variable == 'sint']
      ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
      d_sint = cbind(dNow, setDT(ashObj$result))

      d1 = rbind(d_cost, d_sint)
    }


  # get modified coefficients for k = 1
  if (length(levels(fit$coef_lookup$cond_lev)) > 1){
    d2 = foreach(j = 2:length(levels(fit$coef_lookup$cond_lev)),
                 .combine = rbind) %dopar% {

                   c(ix, iy) %<-% fit$coef_lookup[cond_num %in% c(i,j) & model_comp %like% 'cost']$coef_idx
                   c(jx, jy) %<-% fit$coef_lookup[cond_num %in% c(i,j) & model_comp %like% 'sint']$coef_idx

                   formList = list(~ x1 + x2, ~x3 + x4)
                   idx = c(ix, iy, jx, jy)

                   seMat = getStdErrs(formList, co[, idx], fit$s2.post,
                                      fit$cov.coefficients[idx, idx])

                   statsDT[, cost_k1 := co[, ix] + co[, iy]]
                   statsDT[, sint_k1 := co[jx] + co[, jy]]
                   d3 = melt(statsDT[, .(feature, cost = cost_k1, sint = sint_k1)],
                             measure = c("cost", "sint"),
                             value.name = "estimate")
                   d3[, cond := fit$coef_lookup[cond_num == j]$cond_lev[1L]]
                   d3[, se := as.vector(seMat)]

                   if(isTRUE(combine_dist)){
                     ashObj = d3[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
                     d4 = cbind(d3, setDT(ashObj$result)) } else {

                       dNow = d3[variable == 'cost']
                       ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
                       d_cost = cbind(dNow, setDT(ashObj$result))


                       dNow = d3[variable == 'sint']
                       ashObj = dNow[, ashr::ash(estimate, se, mixcompdist = mixcompdist)]
                       d_sint = cbind(dNow, setDT(ashObj$result))

                       d4 = rbind(d_cost, d_sint)
                     }
                 }


    d = rbind(d1,d2)} else{

      d = d1}

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
                 feature ~ variable+cond, value.var = c("PosteriorMean"))
  coPost = as.matrix(coPost, rownames = 'feature')

  coef_lookup_post = data.table(model_names = colnames(coPost),
                                coef_idx = 1:length(colnames(coPost)))
  coef_lookup_post[,cond := str_extract(model_names, '(?<=_).*')]
  coef_lookup_post[, model_comp := str_extract(model_names, '.*(?=_)')]
  coef_lookup_post[, cond_num := rep(1:length(unique(cond)), 2)]

  c(jx, jy) %<-% coef_lookup_post[cond_num == i]$coef_idx
  dPost = data.table(cond = fit$coef_lookup[cond_num == i]$cond_lev[i],
                     feature = rep(rownames(coRaw), 3),
                     variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coRaw)))
  dPost[variable == 'mesor', estimate_mod := coRaw[,i]]
  dPost[variable == 'amp', estimate_mod := sqrt(coPost[, jx]^2 + coPost[, jy]^2)]
  dPost[variable == 'phase', estimate_mod := atan2(coPost[ jy], coPost[, jx])]

  d0 = merge(dRaw, dPost, by = c("feature", "variable", "cond"))

  #calculate amplitude and phase for k = 1
  if (length(levels(fit$coef_lookup$cond_lev)) >= 2) {
    d1 = foreach(k = 2:length(levels(fit$coef_lookup$cond_lev)),
                 .combine = rbind) %do% {

                   c(kx, ky, kz) %<-% fit$coef_lookup[cond_num == k]$coef_idx
                   dRaw = data.table(cond = fit$coef_lookup[cond_num == k]$cond_lev[i],
                                     feature = rep(rownames(coRaw), 3),
                                     variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coPost)))

                   dRaw[variable == 'mesor',
                        estimate_raw := coRaw[, ix] + coRaw[, kx]]
                   dRaw[variable == 'amp',
                        estimate_raw := sqrt((coRaw[, iy] + coRaw[, ky])^2 + (coRaw[, iz] + coRaw[kz])^2)]
                   dRaw[variable == 'phase',
                        estimate_raw := atan2(coRaw[, iz] + coRaw[, kz], coRaw[, iy] + coRaw[ky])]

                   c(lx, ly) %<-% coef_lookup_post[cond_num == k]$coef_idx

                   dPost = data.table(cond = fit$coef_lookup[cond_num == k]$cond_lev[i],
                                      feature = rep(rownames(coRaw), 3),
                                      variable = rep(c('mesor', 'amp', 'phase'), each = nrow(coRaw)))

                   dPost[variable == 'mesor',
                         estimate_mod := coRaw[, ix] + coRaw[, kx]]
                   dPost[variable == 'amp',
                         estimate_mod := sqrt(coPost[, lx]^2 + coPost[, ly]^2)]
                   dPost[variable == 'phase', estimate_mod := atan2(coPost[, ly], coPost[, lx])]

                   tl = merge(dRaw, dPost, by = c("feature", "variable", "cond"))}

    d = rbind(d0, d1)} else{

      d = d0}

  d[variable == 'phase' & estimate_raw < 0, estimate_raw := estimate_raw +  2*pi]
  d[variable == 'phase' & estimate_raw >= 2*pi, estimate_raw := estimate_raw - 2*pi]
  d[variable == 'phase', estimate_raw := estimate_raw * fit$period / 2 / pi]
  d[variable == 'phase' & estimate_raw > fit$period / 2, estimate_raw := estimate_raw - fit$period/2]

  d[variable == 'phase' & estimate_mod < 0, estimate_mod := estimate_mod + 2*pi]
  d[variable == 'phase' & estimate_mod >=  2*pi, estimate_mod := estimate_mod - 2*pi]
  d[variable == 'phase', estimate_mod := estimate_mod * fit$period / 2 / pi]
  d[variable == 'phase' & estimate_mod > fit$period / 2, estimate_mod := estimate_mod - fit$period/2]

  return(d[])}



getDiffRhythmStats2 =  function(fit, rhythmStats,
                                conds = levels(fit$coef_lookup$cond_lev)[1:2]){

  stopifnot(length(conds) >= 2, all(conds %in% levels(fit$coef_lookup$cond_lev)))

  d1 = rhythmStats[, .(estimate_raw = diff(estimate_raw),
                       estimate_mod = diff(estimate_mod)),
                   by = .(feature, variable)]

  d1[, variable := paste0('diff_', variable)]

  d1[variable == 'diff_phase', estimate_raw := estimate_raw %% fit$period]
  d1[variable == 'diff_phase' & estimate_raw > fit$period / 2,
     estimate_raw := estimate_raw - fit$period]
  d1[variable == 'diff_phase' & estimate_raw < (-fit$period / 2),
     estimate_raw := estimate_raw + fit$period]



  d1[variable == 'diff_phase', estimate_mod := estimate_mod %% fit$period]
  d1[variable == 'diff_phase' & estimate_mod > fit$period / 2,
     estimate_mod := estimate_mod - fit$period]
  d1[variable == 'diff_phase' & estimate_mod < (-fit$period / 2),
     estimate_mod := estimate_mod + fit$period]



  return(d1)}



#' @export

getRhythmTests = function(fit, qvalRhyCutoff = 0.2){
  topTables = getTopTables(fit, qvalRhyCutoff)

  d = merge(topTables$cond_time[, .(feature, pval_diff_rhy = P.Value, qval_diff_rhy = adj.P.Val)],
            topTables$time[, .(feature, pval_rhy = P.Value, qval_rhy = adj.P.Val)],
            by = 'feature', sort = FALSE)


  setorderv(d, c('qval_diff_rhy', 'qval_rhy'))
  return(d[])}




