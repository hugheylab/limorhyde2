# limorhyde2

[![check-deploy](https://github.com/hugheylab/limorhyde2/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/limorhyde2/actions)
[![codecov](https://codecov.io/gh/hugheylab/limorhyde2/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/limorhyde2)
[![Netlify Status](https://api.netlify.com/api/v1/badges/2303634f-911d-4872-85ba-d3a04ed0b952/deploy-status)](https://app.netlify.com/sites/stupefied-engelbart-0482ba/deploys)

`limorhyde2` is an approach to analyze rhythmic, genome-scale data in a way that focuses on effect sizes. For a detailed description of `limorhyde2` along with examples showing its utility on data from one or more conditions, please wait for the preprint.

## Installation

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
