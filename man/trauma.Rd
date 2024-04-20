\name{trauma}
\docType{data}
\alias{trauma}
\title{
Trauma patients from different hospitals
}
\description{
This data set consists of data of 371 trauma patients from three hospitals.
The binary variable \code{mortality} is used as an outcome, and variables \code{age}, \code{sex}, the Injury Severity Score (\code{ISS}, ranging from 1 (low) to 75 (high)) and the Glasgow Coma Scale (\code{GCS}, which expresses the level of consciousness, ranging from 3 (low) to 15 (high)) are used as covariates.
There are three types of hospitals: peripheral hospital without a neuro-surgical unit (\code{Status = 1}), peripheral hospital with a neuro-surgical unit (\code{Status = 2}), and academic medical center (\code{Status = 3}). Originally, the data come from a multi center study collected with a different aim. For educational purposes minor changes have been made, see the references below.
}
\usage{data(trauma)}
\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2024). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>

Draaisma J.M.Th, de Haan A.F.J., Goris R.J.A. (1989). \emph{Preventable Trauma Deaths in the Netherlands - A prospective Multicentre Study}, {The journal of Trauma}, Vol. 29(11), 1552-1557.
}
\keyword{datasets}
