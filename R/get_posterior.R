#' @export
getPosteriorSamples = function(fit, nPosteriorSamples = 100, overwrite = FALSE) {

  stopifnot(!is.null(fit$mashFit),
            isTRUE(overwrite) || is.null(fit$mashPosteriorSamples),
            length(nPosteriorSamples) == 1L,
            is.numeric(nPosteriorSamples),
            nPosteriorSamples >= 10)

  nPostSamps = round(nPosteriorSamples)

  mp = mashr::mash_compute_posterior_matrices(
    fit$mashFit, fit$mashData, algorithm.version = 'R',
    posterior_samples = nPostSamps)

  co = fit$coefficients
  coArray = array(rep(co, nPostSamps), dim = c(dim(co), nPostSamps),
                  dimnames = dimnames(co))
  coArray[, fit$mashIdx, ] = mp$PosteriorSamples

  fit$mashPosteriorSamples = coArray
  return(fit)}


# TODO
# #' @export
# getPosteriorUncertainty = function(posteriorStats, ) {
#   0
#   return(d)} # data.table with one row per feature per statistic
