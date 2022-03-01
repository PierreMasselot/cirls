# cirls: Constrained iteratively reweighted least-squares

The package `cirls` provides routines to fit Generalized Linear Models (GLM) with coefficients subject to linear constraints, through a constrained iteratively reweighted least-squares algorithm. 

## Installation

The package `cirls` is currently *under devopment* and not yet submitted to CRAN. The current version can be installed through the `devtools` package as

```R
install_github("PierreMasselot/cirls")
```

and can then be loaded as usual `library(cirls)`.

## Usage

The central function of the package is `cirls.fit` meant to be passed through the `method` argument of the `glm` function. The user is also expected to pass a either constraint matrix or a list of constraint matrices through the `Cmat` argument. Typical usage is then:

```R
consres <- glm(y ~ x + z, method = "cirls.fit", Cmat = Cmat)
```

or equivalently

```R
consres <- glm(y ~ x + z, method = "cirls.fit", Cmat = list(x = cx, z = cz))
```

## References

*To come*
