## R CMD check results

### Local check
`devtools::check()` result:

  0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### Online check
`devtools::check_rhub()` result:

  > checking CRAN incoming feasibility ... NOTE
    New submission
    Possibly misspelled words in DESCRIPTION:
      Rhythmicity (3:17, 3:46)
    
    
    Maintainer: 'Jake Hughey <jakejhughey@gmail.com>'
  
  > checking examples ... NOTE
    Examples with CPU (user + system) or elapsed time > 5s
                        user system elapsed
    getDiffRhythmStats  6.33   0.22    6.67
    getStatsIntervals   5.04   0.11    5.21
    getPosteriorSamples 6.78   0.09    7.40
    getRhythmStats      4.28   0.10    5.56
  
  > checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
      'lastMiKTeXException'
  
  0 errors ✓ | 0 warnings ✓ | 3 notes x

Notes:

  - This is the first time this package is being submitted to CRAN and the word is not misspelled, so ignore this note.
  
  - The examples here are taking a little over 5 seconds for the functions listed, but are still under 10 seconds.
  
  - The third note only occurs on the Windows Server rhub environment, and from what I have seen about these types of notes they do not occur when building and checking on CRAN.

You can also see the results of R CMD check on Windows, Linux, and MacOS by going to the GitHub Actions run labeled `check-deploy` [here](https://github.com/hugheylab/limorhyde2/actions).

## Resubmission notes
The initial submission was automatically rejected due to the examples taking over 10 seconds to run. To address this we did the following

  - We reduced the data used in the examples to reduce processing time.
  - Examples are now under 10 seconds.

## Downstream dependencies
There are no downstream dependencies for limorhyde2.
