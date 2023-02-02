## R CMD check results

### Local

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### R-hub (Windows)

  Found the following files/directories:
    'lastMiKTeXException'

  0 errors ✔ | 0 warnings ✔ | 1 note ✖

See [here](https://builder.r-hub.io/status/limorhyde2_0.1.0.tar.gz-ea0800b3d4064efe97b2b5e26515847e).

### GitHub Actions (Mac, Ubuntu, Windows)

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

See [here](https://github.com/hugheylab/limorhyde2/actions/runs/4018747414).

## Changes from current CRAN release

* Changed default for credible interval mass to 0.9.
* Fixed `getExpectedMeas()` for single condition and no covariates.
* Added `isAlreadyInParallel()` function that checks if already running in parallel.
* Added check to `getRhythmStats()` to prevent running in parallel if already in parallel.
* Fixed `mergeMeasMeta()` for RNA-seq data.
* Updated code style.
* Enabled spline fit with fewer knots and customizable degree.
