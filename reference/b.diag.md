# Create a Block Diagonal Matrix

Construct a block diagonal matrix using multiple given block matrices.

## Usage

``` r
b.diag(..., fill = 0)
```

## Arguments

- ...:

  individual matrices or one list of matrices.

- fill:

  non-block-diagonal elements. Default is \\0\\.

## Details

Avoid combining matrices and lists for the `...` argument.

`b.diag()` covers the arguments of type "character".

If a *sparse* matrix needed, run the following:

[`library(Matrix); Matrix(b_diag, sparse = TRUE)`](https://Matrix.R-forge.R-project.org)

where `b_diag` is the matrix returned by `b.diag()`.

## Value

`b.diag()` returns a block diagonal matrix obtained by combining the
arguments.

## Author

Hassan Pazira  
Maintainer: Hassan Pazira <h.pazira@arq.org>

## Examples

``` r
b.diag(1, matrix(1:3, 3,4), diag(3:2))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    1    0    0    0    0    0    0
#> [2,]    0    1    1    1    1    0    0
#> [3,]    0    2    2    2    2    0    0
#> [4,]    0    3    3    3    3    0    0
#> [5,]    0    0    0    0    0    3    0
#> [6,]    0    0    0    0    0    0    2

b.diag(matrix(1:6, 2), as.character(2))
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    3    5    0
#> [2,]    2    4    6    0
#> [3,]    0    0    0    2

lists <- list(1, 2:3, diag(4:6), 7, cbind(8,9:12), 13:15)
b.diag(lists)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]    1    0    0    0    0    0    0    0    0
#>  [2,]    0    2    0    0    0    0    0    0    0
#>  [3,]    0    3    0    0    0    0    0    0    0
#>  [4,]    0    0    4    0    0    0    0    0    0
#>  [5,]    0    0    0    5    0    0    0    0    0
#>  [6,]    0    0    0    0    6    0    0    0    0
#>  [7,]    0    0    0    0    0    7    0    0    0
#>  [8,]    0    0    0    0    0    0    8    9    0
#>  [9,]    0    0    0    0    0    0    8   10    0
#> [10,]    0    0    0    0    0    0    8   11    0
#> [11,]    0    0    0    0    0    0    8   12    0
#> [12,]    0    0    0    0    0    0    0    0   13
#> [13,]    0    0    0    0    0    0    0    0   14
#> [14,]    0    0    0    0    0    0    0    0   15
identical(b.diag(lists), b.diag(lapply(lists, as.matrix)))
#> [1] TRUE

b.diag(replicate(3, matrix(round(rnorm(9)), 3, 3), simplify=FALSE))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#>  [1,]   -1    0   -2    0    0    0    0    0    0
#>  [2,]    0    1    0    0    0    0    0    0    0
#>  [3,]   -2    1    0    0    0    0    0    0    0
#>  [4,]    0    0    0    0    2   -2    0    0    0
#>  [5,]    0    0    0   -1   -2   -1    0    0    0
#>  [6,]    0    0    0    1    1    0    0    0    0
#>  [7,]    0    0    0    0    0    0    1    0    2
#>  [8,]    0    0    0    0    0    0   -1   -1    0
#>  [9,]    0    0    0    0    0    0    0    1   -1
```
