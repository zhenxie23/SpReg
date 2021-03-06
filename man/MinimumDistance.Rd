% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MinimumDistance.R
\name{MinimumDistance}
\alias{MinimumDistance}
\alias{summary.MinimumDistance}
\title{Estimation and Inference for Minimum-Distance Estimator}
\usage{
MinimumDistance(DatR, VarR, DatY, VarY, variogram.model,
                       MD.starting.value, cutoff.R, cutoff.YR,
                       start.value.method = 2, projected = FALSE)

\method{summary}{MinimumDistance}(object, ...)
}
\arguments{
\item{DatR}{explanatory variable R, a spatial object, see \link[sp:coordinates]{coordintes()}}

\item{VarR}{name of variable R}

\item{DatY}{outcome variable Y, a spatial object, see \link[sp:coordinates]{coordintes()}}

\item{VarY}{name of variable Y}

\item{variogram.model}{variogram model type, e.g. "Exp", "Sph", "Gau", "Mat"}

\item{MD.starting.value}{the starting point of parameters}

\item{cutoff.R}{cutoff for sample variogram of variable R}

\item{cutoff.YR}{cutoff for sample cross variogram of variable R and Y}

\item{start.value.method}{fitting method, see \link[gstat:fit.variogram]{fit.variogram()}}

\item{projected}{logical; if FALSE, data are assumed to be unprojected, meaning decimal longitude/latitude. For projected data, Euclidian distances are computed, for unprojected great circle distances(km) are computed.}

\item{object}{class \code{MinimumDistance} objects.}
}
\value{
\item{\code{num.obs}}{the number of observations}
        \item{\code{vario.par.point.est}}{point estimates for variogram parameters(psill, range)}
        \item{\code{vario.par.var.mat}}{estimated variance-covariance matrix for variogram parameters(psill, range)}
        \item{\code{md.point.est}}{point estimates for Min-Dist estimator}
        \item{\code{md.var.mat}}{estimated aymptotic variance for Min-Dist estimator}
}
\description{
This function implements the esmation and large sample inference for
             Minimum-Distance Estimator.
}
\section{Methods (by generic)}{
\itemize{
\item \code{summary}: \code{summary} method for class "\code{MinimumDistance}".
}}

\seealso{
\link{sp}, \link{gstat}
}
