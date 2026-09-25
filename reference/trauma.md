# Trauma patients from different hospitals

This data set consists of data of 371 trauma patients from three
hospitals. The binary variable `mortality` is used as an outcome, and
variables `age`, `sex`, the Injury Severity Score (`ISS`, ranging from 1
(low) to 75 (high)) and the Glasgow Coma Scale (`GCS`, which expresses
the level of consciousness, ranging from 3 (low) to 15 (high)) are used
as covariates. There are three types of hospitals: peripheral hospital
without a neuro-surgical unit (`Status = 1`), peripheral hospital with a
neuro-surgical unit (`Status = 2`), and academic medical center
(`Status = 3`). Originally, the data come from a multi center study
collected with a different aim. For educational purposes minor changes
have been made, see the references below.

## Usage

``` r
data(trauma)
```

## References

Jonker M.A., Pazira H. and Coolen A.C.C. (2024). *Bayesian federated
inference for estimating statistical models based on non-shared
multicenter data sets*, *Statistics in Medicine*, 43(12): 2421-2438.
\<https://doi.org/10.1002/sim.10072\>

Draaisma J.M.Th, de Haan A.F.J., Goris R.J.A. (1989). *Preventable
Trauma Deaths in the Netherlands - A prospective Multicentre Study*, The
journal of Trauma, Vol. 29(11), 1552-1557.
