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
