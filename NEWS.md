# limorhyde2 0.0.10
* Fixed `mergeMeasMeta()` for RNA-seq data.

# limorhyde2 0.0.9
* Enabled spline fit with fewer knots and customizable degree.

# limorhyde2 0.0.8
* Matched style to lab standard.

# limorhyde2 0.0.7
* Made data smaller and examples faster.

# limorhyde2 0.0.6
* Switched to using updated `ashr` from CRAN.

# limorhyde2 0.0.4
* Added more stringent checking for correspondence between measurements matrix and metadata.

# limorhyde2 0.0.3
* Made `mash` silent.
* Added checking for parallel backend, to avoid warning messages.

# limorhyde2 0.0.2
* Moved data to data/, to make examples and vignettes simpler.
* Added examples.

# limorhyde2 0.0.1
* Updated vignette for analyzing RNA-seq data.

# limorhyde2 0.0.0.9033
* Renamed `mean_value` to `mesor` and added differential rhythm statistics for mean mesor and mean amp.

# limorhyde2 0.0.0.9032
* Added differential rhythm statistic `diff_rhy_dist`.

# limorhyde2 0.0.0.9031
* Cleaned up calculation of credible intervals for rhythm amplitudes.

# limorhyde2 0.0.0.9030
* Updated model parameterization and mashing to treat all conditions fairly.

# limorhyde2 0.0.0.9029
* Updated calculation of credible intervals for rhythm amplitudes.

# limorhyde2 0.0.0.9028
* Added vignette for RNA-seq data.

# limorhyde2 0.0.0.9027
* Added ability to fit models using `DESeq2`.

# limorhyde2 0.0.0.9026
* Made calculation of rms-based stats optional.
* Simplified tests to use snapshotting.

# limorhyde2 0.0.0.9025
* Changed `assertLogical()` to `assertFlag()` for checking inputs.

# limorhyde2 0.0.0.9024
* Fixed vignette building on Windows machines

# limorhyde2 0.0.0.9023
* Fixed factor weirdness between R versions

# limorhyde2 0.0.0.9022
* Output for `getDiffRhythmStats()` now includes all pairs of conditions by default

# limorhyde2 0.0.0.9021
* Fixed bug so getDiffRhythmStats uses all shifted models

# limorhyde2 0.0.0.9020
* Replaced `globalVariables()` with assigning variables to `NULL` in each function

# limorhyde2 0.0.0.9019
* Clarified wording in vignettes

# limorhyde2 0.0.0.9018
* Added function to merge measurements and metadata
* Added vignettes

# limorhyde2 0.0.0.9017
* Fixed bug in getStatsIntervals

# limorhyde2 0.0.0.9016
* Cleaned up logic around nKnots for cosinor
* Fixed setting attributes for data.tables
* Fixed input check in getDiffRhythmStats

# limorhyde2 0.0.0.9015
* Revamped input checking to use the checkmate package
