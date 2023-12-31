% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rdirichlet_via_gamma}
\alias{rdirichlet_via_gamma}
\title{Random sample from a Dirichlet distribution through gamma}
\usage{
rdirichlet_via_gamma(n, alpha)
}
\arguments{
\item{n}{integer : number of items to be sampled}

\item{alpha}{scalar, with \eqn{\alpha > 0} : concentration parameter of the
distribution}
}
\value{
A numeric matrix of shape \eqn{n \times \mathrm{length}(\alpha)}
storing an $n$-sample drawn from a Dirichlet distribution with parameters
\eqn{\alpha}.
}
\description{
Random sample from a Dirichlet distribution through gamma
}
