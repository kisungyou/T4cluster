% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_gen3S.R
\name{gen3S}
\alias{gen3S}
\title{Generate from Three 5-dimensional Subspaces in 200-dimensional space.}
\usage{
gen3S(n = 50, var = 0.3)
}
\arguments{
\item{n}{the number of data points sampled from each subspace (default: 50).}

\item{var}{degree of Gaussian noise (default: 0.3).}
}
\value{
a named list containing with :\describe{
\item{data}{an \eqn{(3*n\times 3)} data matrix.}
\item{class}{length-\eqn{3*n} vector for class label.}
}
}
\description{
Generate from Three 5-dimensional Subspaces in 200-dimensional space.
}
\examples{
\donttest{
## a toy example
tester = gen3S(n=100)
data   = tester$data
label  = tester$class
}

}
\references{
\insertRef{wang_efficient_2011}{T4cluster}
}
\concept{data}
