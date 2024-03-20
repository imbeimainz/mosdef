# mosdef

The goal of `mosdef` is to provides functionality to run a number of tasks in the differential expression analysis workflow. 

This encompasses the most widely used steps, from running various enrichment analysis tools with a unified interface to creating plots and beautifying table components linking to external websites and databases. 

This streamlines the generation of comprehensive analysis reports.

`mosdef` can be found on Bioconductor
(<https://www.bioconductor.org/packages/mosdef>).

## Installation

You can install the development version of `mosdef` from GitHub with:

``` r
library("remotes")
remotes::install_github("imbeimainz/mosdef", 
                        dependencies = TRUE, build_vignettes = TRUE)
```

## Usage overview

You can find the rendered version of the documentation of `mosdef` at
the project website <https://imbeimainz.github.io/mosdef>, created with `pkgdown`.

## Development

If you encounter a bug, have usage questions, or want to share ideas and
functionality to make this package better, feel free to file an
[issue](https://github.com/imbeimainz/mosdef/issues).

## Code of Conduct

Please note that the mosdef project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.

## License

MIT

