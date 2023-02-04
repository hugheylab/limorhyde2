# limorhyde2

[![check-deploy](https://github.com/hugheylab/limorhyde2/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/limorhyde2/actions)
[![codecov](https://codecov.io/gh/hugheylab/limorhyde2/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/limorhyde2)
[![Netlify Status](https://api.netlify.com/api/v1/badges/2303634f-911d-4872-85ba-d3a04ed0b952/deploy-status)](https://app.netlify.com/sites/stupefied-engelbart-0482ba/deploys)
[![CRAN Status](https://www.r-pkg.org/badges/version/limorhyde2)](https://cran.r-project.org/package=limorhyde2)
[![drat version](https://raw.githubusercontent.com/hugheylab/drat/gh-pages/badges/limorhyde2_drat_badge.svg)](https://github.com/hugheylab/drat/tree/gh-pages/src/contrib)

`limorhyde2` is an approach to analyze rhythmic, genome-scale data in a way that focuses on effect sizes. For a detailed description of `limorhyde2` along with examples showing its utility for transcriptome data, please check out the [preprint](https://doi.org/10.1101/2023.02.02.526897).

## Installation

### Option 1: CRAN

```r
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'limma'))
install.packages('limorhyde2')
```

### Option 2: Hughey Lab Drat Repository

1. Install [`BiocManager`](https://cran.r-project.org/package=BiocManager).

    ```r
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    ```

1. If you use RStudio, go to Tools → Global Options... → Packages → Add... (under Secondary repositories), then enter:

    - Name: hugheylab
    - Url: https://hugheylab.github.io/drat/

    You only have to do this once. Then you can install or update the package by entering:

    ```r
    BiocManager::install('limorhyde2')
    ```

    Alternatively, you can install or update the package by entering:

    ```r
    BiocManager::install('limorhyde2', site_repository = 'https://hugheylab.github.io/drat/')
    ```

## Usage

See the examples and detailed guidance in the [reference documentation](https://limorhyde2.hugheylab.org/reference/index.html).
