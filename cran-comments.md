## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

  > checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Jake Hughey <jakejhughey@gmail.com>'

    New submission

    Possibly misspelled words in DESCRIPTION:
      Rhythmicity (3:17, 3:46)
      mutivariate (11:22)

  > checking examples ... NOTE
    Examples with CPU (user + system) or elapsed time > 5s
                        user system elapsed
    getStatsIntervals   7.47   0.17    7.64
    getPosteriorSamples 7.24   0.18    7.40

  > checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
      'lastMiKTeXException'

  0 errors ✓ | 0 warnings ✓ | 3 notes x

Notes:

  - This is the first time this package is being submitted to CRAN and the words are not misspelled, so ignore this note.
  
  - The examples here are taking a little over 7 seconds for the two functions listed. If this is a problem, let us know.
  
  - The third note only occurs on the Windows Server rhub environment, and from what I have seen about these types of notes they do not occur when building and checking on CRAN.

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/limorhyde2/actions).

## Downstream dependencies
There are no downstream dependencies for limorhyde2.
