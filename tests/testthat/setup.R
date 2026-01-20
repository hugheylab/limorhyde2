library('glue')

foreach::registerDoSEQ()
dataDir = 'data'

snapshot = function(xObs, path) {
  if (file.exists(path)) {
    # xExp = qs::qread(path)
    xExp = readRDS(path)
  } else {
    # qs::qsave(xObs, path)
    saveRDS(xObs, path)
    xExp = xObs}
  return(xExp)}
